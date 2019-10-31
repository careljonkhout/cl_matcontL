function [V, reduced_jacobian, delta_q_gamma, delta_q_r, G_delta_q_r, ...
          phases_0, phases_T_i, period, active_par_val] = ...
            compute_reduced_jacobian(x)

  global cds contopts
  [phases_0, period, parameters] = ...
    NewtonPicard.MultipleShooting.extract_phases_period_and_parameters(x);
  cds.phases_0 = phases_0;
  active_par_val = parameters{cds.ActiveParams};
  
  if ~ cds.using_cvode
    % Since the result of the next time-integration will be reused many times,
    % we set the tolerances a bit tighter than the rest.
    integration_opt = odeset(...
      'AbsTol',      contopts.integration_abs_tol / 10,    ...
      'RelTol',      contopts.integration_rel_tol / 10    ...
    );

    if ~ isempty(cds.jacobian_ode)
      integration_opt = odeset(integration_opt, ...
      'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}));
    end
  end


  m          = cds.nMeshIntervals;
  phases_T_i = zeros(cds.nphases, m);
  
  integrator            = cds.integrator;
  dydt_ode              = @(t,y) cds.dydt_ode(t, y, parameters{:});
  mesh                  = cds.mesh;
  % compute length of each mesh interval:
  delta_t               = period * diff(mesh);
  
  if cds.using_cvode
    [~,y] = feval(integrator, ...
      'initial_point',   phases_0(:,1), ...
      't_values',        [0 delta_t(1)], ...
      'ode_parameters',  cell2mat(parameters), ...
      'abs_tol',         contopts.integration_abs_tol, ...
      'rel_tol',         contopts.integration_rel_tol);
    phases_T_i(:,1) = y(end,:)';
  else
    cds.orbits(1) = feval(integrator,...
            @(t, y) feval(dydt_ode, t, y), ...
            [0 period], ...
            phases_0(:,1), integration_opt);
    phases_T_i(:,1) = deval(cds.orbits(1), delta_t(1));
  end
  

  if contopts.contL_ParallelComputing && cds.using_cvode
    error('parallel computing with cvode not yet implemented for multiple shooting')
  end
  parfor_failed = false;
  try 
    if contopts.contL_ParallelComputing
      parfor i=2:m
        % The function monodromy_map cannot be used here, since it depends on
        % the global variable cds, and global variables are not copied so the
        % the workspace of the workers that parfor uses.
        orbits(i) = feval(integrator,...
          @(t, y) feval(dydt_ode, t, y), ...
          [0 delta_t(i)], ...
          phases_0(:,i), integration_opt);
        phases_T_i(:,i) = deval(orbits(i), delta_t(i));
      end
      cds.orbits(2:m) = orbits(2:m);
    end
    catch my_error
    if (strcmp(my_error.identifier, ...
                                'MATLAB:remoteparfor:AllParforWorkersAborted'))
      % Something went wrong with the parfor workers.
      % We try again with ordinary for.
      print_diag(0, 'Parfor aborted, retrying with ordinary for.\n');
      parfor_failed = true;
    else
      % in case of some other error, we want to know about it, so we rethrow the
      % error
      rethrow(my_error)
    end
  end
  
  if ~ contopts.contL_ParallelComputing || parfor_failed
   
    for i=2:m
      print_diag(6,'orbit %d\n',i)
      if cds.using_cvode
        [~, y] = feval(integrator, ...
          'initial_point',   phases_0(:,i), ...
          't_values',        [0 delta_t(i)], ...
          'ode_parameters',  cell2mat(parameters), ...
          'abs_tol',         contopts.integration_abs_tol, ...
          'rel_tol',         contopts.integration_rel_tol);
        phases_T_i(:,i) = y(end,:)';
      else
        cds.orbits(i) = feval(integrator,...
          @(t, y) feval(dydt_ode,t, y), ...
          [0 delta_t(i)], ...
          phases_0(:,i), integration_opt);

        phases_T_i(:,i) = deval(cds.orbits(i), delta_t(i));
      end
    end
  end

  V    = cell(m,1);
  V{1} = compute_subspace(1, period, parameters);
  V{2} = compute_subspace(2, period, parameters);
  MV   = cell(m,1);
  M    = @NewtonPicard.MultipleShooting.monodromy_map; 
  
  for j = 1 : size(V{1}, 2)
    MV{1}(:, j) = M(i, V{1}(:, j), delta_t(1), parameters);
  end
  
  for i = 2 : m % m == cds.nMeshIntervals
    print_diag(6,'G %d\n',i);
    for j = 1 : size(V{i}, 2)
      MV{i}(:, j) = M(i, V{i}(:, j), delta_t(i), parameters);
    end
    if i < m
      V_i_plus_one = orth(MV{i});
      if size(V_i_plus_one, 2) >= cds.preferred_basis_size
        V{i+1} = V_i_plus_one;
      else
        print_diag(1, ...
                 'subspace is rank deficient. computing subspace using eigs\n');
        V{i + 1} = compute_subspace(i + 1, period, parameters);
      end
    end
  end
  
  cds.V = V;
  
  % we compute the number of nonzero's (nnz):
  block_size = cds.preferred_basis_size + 1 ;
  nnz = m * block_size^2;      % blocks on diagonal;
  nnz = nnz + m * block_size;  % identity matrices on super-block-diagonal
  reduced_jac_size = 0;
  for i = 1 : m
    reduced_jac_size = reduced_jac_size + size(V{i}, 2);
  end
  cds.reduced_jac_size = reduced_jac_size;
  
  % the variable name F_0_pp corresponds with L^0_pp in Lust
  F_0_pp = spalloc(reduced_jac_size, reduced_jac_size, nnz);
  
  row_offset = 0;
  col_offset = 0;
  for i=1:m % note that: m == cds.nMeshIntervals
    ni = next_index_in_cycle(i,m);
    row_indices = row_offset + (1 : size(V{ni}, 2));
    col_indices = col_offset + (1 : size(V{ i}, 2));
    F_0_pp(row_indices, col_indices) = V{ni}' * MV{i}; %#ok<SPRIX>
    % Indexation might be slow, but it the most time-consuming part of the
    % algorithm. (The most time consuming part of the algorithm is
    % time-integration). Therefore we ignore the <SPRIX> code analyzer warnings.
    col_offset = col_offset + size(V{ i}, 2);
    row_offset = row_offset + size(V{ni}, 2);
  end

  % add the negative identity matrices to F_0_pp
  row_offset = 0;
  col_offset = size(V{1}, 2);
  for i = 1 : m-1
    ni          = next_index_in_cycle(i, m);
    eye_size    = size(V{ni}, 2);
    row_indices = row_offset + (1 : eye_size);
    col_indices = col_offset + (1 : eye_size);
    F_0_pp(row_indices, col_indices) = - eye(eye_size); %#ok<SPRIX>
    row_offset  = row_offset + eye_size;
    col_offset  = col_offset + eye_size;
  end
  eye_size = size(V{m}, 2);
  F_0_pp( end-eye_size+1 : end,  1 : eye_size) = - eye(eye_size); 
  
  
  V_T__b_T = zeros(reduced_jac_size, 1); % V_p^T b_T      in \cite{lust-phd}
  
  row_offset = 0;
  for i=1:m
    ni                = next_index_in_cycle(i,m);
    indices           = row_offset + (1 : size(V{ni},2));
    b_T_i             = cds.dydt_ode(0, phases_T_i(:,i), parameters{:}) * ...
                                delta_t(i) / period;
    V_T__b_T(indices) = V{ni}' * b_T_i;
    row_offset        = row_offset + size(V{ni},2);
  end
  

  rhs_delta_q_gamma = zeros(cds.nphases,m);
 
  V_T__b_g = zeros(reduced_jac_size, 1);
  
  row_offset = 0;
  for i=1:m
    b_g_i                  = NewtonPicard.compute_d_phi_d_p(phases_0(:,i), ...
                                           delta_t(i), parameters);
    ni                     = next_index_in_cycle(i,m);
    indices                = row_offset + (1 : size(V{ni},2));
    V_T__b_g(indices)      = V{ni}' * b_g_i;
    
    rhs_delta_q_gamma(:,i) = b_g_i - V{ni} * V_T__b_g(indices);
    
    row_offset             = row_offset + size(V{ni},2);
  end
  
  
  % the r in q_r means residual
  
  rhs_delta_q_r = zeros(cds.nphases,m);
  
  
  for i=1:m % m == cds.nMeshIntervals
    ni                 = next_index_in_cycle(i,m);
    r                  = phases_T_i(:,i) - phases_0(:,ni);
    rhs_delta_q_r(:,i) = r - V{ni} * (V{ni}' * r);
  end
  
  solve_Q_system = @NewtonPicard.MultipleShooting.solve_Q_system;
  
  [delta_q_gamma, G_delta_q_gamma] = solve_Q_system(V, rhs_delta_q_gamma, ...
                                            delta_t, parameters);

  [delta_q_r    , G_delta_q_r    ] = solve_Q_system(V, rhs_delta_q_r, ...
                                            delta_t, parameters);
 
  jac_1_3 = V_T__b_g;
  
  row_offset = 0;
  for i=1:m % m == cds.nMeshIntervals
    ni                = next_index_in_cycle(i,m);
    indices           = row_offset + (1 : size(V{ni},2));
    jac_1_3(indices)  = jac_1_3(indices) + V{ni}' * G_delta_q_gamma(:,i);
    row_offset        = row_offset + size(V{ni},2);
  end

  d_s_d_T      = 0;
  d_s_d_gamma  = 0;
  
  
  jac_2_1                         = zeros(1, reduced_jac_size);
  jac_2_1( 1, 1 : size(V{1}, 2) ) = cds.previous_dydt_0' * V{1};
  
  jac_2_3 = d_s_d_gamma + cds.previous_dydt_0' * delta_q_gamma(:,1);
 
  reduced_jacobian = [
    F_0_pp     V_T__b_T     jac_1_3;
    jac_2_1    d_s_d_T      jac_2_3;
  ];
  
end

  
function V = compute_subspace(i, period, parameters)
  global cds
  
  n_eigenvalues = min(cds.nphases,cds.p + 1);
  
  if i == 1
    n_eigenvalues = min(cds.nphases, n_eigenvalues + 2);
  end
  
  [eigenvectors, eigenvalues, no_convergence] = eigs( ...
    @(x) NewtonPicard.MultipleShooting.monodromy_map( ...
       i, x, period, parameters), ...
    cds.nphases, ...
    n_eigenvalues);


  if no_convergence
    V = [];
    fprintf(['Newton_Picard_Correction.m:', ...
      ' eigenvalues of monodromy matrix did not converge.\n'])
    return
  end

  eigenvalues = diag(eigenvalues);
  basis = zeros(cds.nphases, cds.p + 1);

  i = 1;

  while i <= n_eigenvalues
    basis(:, i) = real(eigenvectors(:,i));
    i = i + 1;
    if max(abs(imag(eigenvalues(i-1)))) > 1e-14
      basis(:, i) = imag(eigenvectors(:,i-1));
      i = i + 1;
    end
  end
  V = orth(basis(:,1:i-1));
end


function index = next_index_in_cycle(i,m)
  index = i + 1;
  if index == m + 1
    index = 1;
  end
end
