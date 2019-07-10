function [V, reduced_jacobian, delta_q_gamma, delta_q_r, G_delta_q_r, ...
          phases_0, phases_T_i, period, active_par_val] = ...
            compute_reduced_jacobian(x)

  global cds contopts
  [phases_0, period, parameters] = ...
    NewtonPicard.MultipleShooting.extract_phases_period_and_parameters(x);
  
  active_par_val = parameters{cds.ActiveParams};
  
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


  m          = cds.nMeshIntervals;
  phases_T_i = zeros(cds.nphases,m);
  
  integrator            = cds.integrator;
  dydt_ode              = @(t,y) cds.dydt_ode(t, y, parameters{:});
  mesh                  = cds.mesh;
  % compute length of each mesh interval:
  delta_t               = period * diff(mesh);
  
  
  cds.orbits(1) = feval(integrator,...
          @(t, y) feval(dydt_ode, t, y), ...
          [0 period], ...
          phases_0(:,1), integration_opt);
  phases_T_i(:,1) = deval(cds.orbits(1), delta_t(1));
  
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
  catch error
    if (strcmp(error.identifier,'MATLAB:remoteparfor:AllParforWorkersAborted'))
      % Something went wrong with the parfor workers.
      % We try again with ordinary for.
      print_diag(0, 'Parfor aborted, retrying with ordinary for.\n');
      parfor_failed = true;
    else
      % in case of some other error, we want to know about it
      rethrow(error)
    end
  end
  
  if ~ contopts.contL_ParallelComputing || parfor_failed
    for i=2:m
      print_diag(6,'orbit %d\n',i)
      % The function monodromy_map cannot be used here, since it depends on
      % the global variable cds, and global variables are not copied so the
      % the workspace of the workers that parfor uses.
      cds.orbits(i) = feval(integrator,...
        @(t, y) feval(dydt_ode,t, y), ...
        [0 delta_t(i)], ...
        phases_0(:,i), integration_opt);

      phases_T_i(:,i) = deval(cds.orbits(i), delta_t(i));
    end
  end
  
  
  
  cds.p = cds.preferred_basis_size;

  if ~ isfield(cds, 'V') || true
    V1             = compute_subspace(1, period, parameters);
    basis_size     = size(V1,2);
    cds.basis_size = basis_size;
    V              = zeros(cds.nphases, basis_size, m);
    V(:,:,1)       = V1;
    for i=2:m % m == cds.nMeshIntervals
      print_diag(6,'G %d\n',i);
      for j = 1:size(V,2)
        V(:,j,i) = NewtonPicard.MultipleShooting.monodromy_map(...
          i-1, V(:,j,i-1), delta_t(i-1), parameters);
      end
      V(:,:,i) = orth(V(:,:,i));
    end
  else
    V = NewtonPicard.MultipleShooting.continue_subspaces(delta_t, parameters);
  end
  cds.V = V;
  basis_size     = size(V,2);
  cds.basis_size = basis_size;
  

  MV = zeros(cds.nphases,basis_size,cds.nMeshIntervals);
 
  int_opt = odeset(...
    'AbsTol', contopts.integration_abs_tol,    ...
    'RelTol', contopts.integration_rel_tol );

  % the variable name F_0_pp corresponds with L^0_pp in Lust
  % we compute the number of nonzero's (nnz):
  nnz = m * basis_size * basis_size; % blocks on diagonal;
  nnz = nnz + m * basis_size;       % identity matrices on superdiagonal
  F_0_pp = spalloc(m*basis_size,m*basis_size,nnz);
  
  for i=1:m % note that: m == cds.nMeshIntervals
%     int_opt = odeset(int_opt, ...
%       'Jacobian', @(t,y) feval(cds.jacobian_ode, ...
%                     t,   deval(cds.orbits(i),t), parameters{:}));
%     dydt_partial_mon = @(t, y) ...
%       cds.jacobian_ode(t,deval(cds.orbits(i),t),parameters{:}) * y;
    

    for j=1:basis_size
      %[~,orbit] = ...
      %  cds.integrator(dydt_partial_mon, [0 delta_t(i)], V(:,j,i), int_opt);
      %MV(:,j,i) = orbit(end,:)';
      MV(:,j,i) = NewtonPicard.MultipleShooting.monodromy_map( ...
                                          i, V(:,j,i), delta_t(i), parameters);
    end

    indices = basis_size*(i-1) + (1:basis_size);
    ni = next_index_in_cycle(i,m);
    F_0_pp(indices,indices)= V(:,:,ni)' * MV(:,:,i); %#ok<SPRIX>
    % Indexation might be slow, but it is not the bottleneck in the code. (The
    % bottle neck is time-integration). Therefore we ignore the <SPRIX> code
    % analyzer warnings.
    
    if i < m
      F_0_pp(indices, indices + basis_size) = - eye(basis_size); %#ok<SPRIX>
    else
      F_0_pp(indices, 1:basis_size)         = - eye(basis_size); %#ok<SPRIX>
    end
  end
  
  
  
  V_T__b_T = zeros(basis_size * m,1); % V_p^T b_T      in \cite{lust-phd}
  
  indices = 1:basis_size;
  for i=1:m
    b_T_i             = cds.dydt_ode(0, phases_T_i(:,i), parameters{:}) * ...
                          delta_t(i) / period;
    V_T__b_T(indices) = V(:,:,next_index_in_cycle(i,m))' * b_T_i;
    indices           = indices + basis_size;
  end
  

  rhs_delta_q_gamma = zeros(cds.nphases,m);
 
  V_T__b_g = zeros(basis_size * m,1);
 
  
  indices = 1:basis_size;
  for i=1:m
    b_g_i                  = NewtonPicard.compute_d_phi_d_p(phases_0(:,i), ...
                                           delta_t(i), parameters);
    
    ni                     = next_index_in_cycle(i,m);
    V_T__b_g(indices)      = V(:,:,ni)' * b_g_i;
    
    rhs_delta_q_gamma(:,i) = b_g_i - V(:,:,ni) * V_T__b_g(indices);
    indices                = indices + basis_size;
  end
  
  
  % the r in q_r means residual
  
  rhs_delta_q_r = zeros(cds.nphases,m);
  
  for i=1:m % m == cds.nMeshIntervals
    ni                 = next_index_in_cycle(i,m);
    r                  = phases_T_i(:,i) - phases_0(:,ni);
    rhs_delta_q_r(:,i) = r - V(:,:,ni) * V(:,:,ni)' * r;
  end
  
  [delta_q_gamma, G_delta_q_gamma] = ...
    NewtonPicard.MultipleShooting.solve_Q_system(V, ...
    rhs_delta_q_gamma, delta_t, parameters);

  [delta_q_r    , G_delta_q_r]     = ...
    NewtonPicard.MultipleShooting.solve_Q_system(V, ...
    rhs_delta_q_r, delta_t, parameters);
  
  
  V_T_d_phi_d_T = zeros(basis_size*m,1);
  for i=1:m % m == cds.nMeshIntervals
    indices = (i-1)*basis_size + (1:basis_size);
    ni = next_index_in_cycle(i,m);
    V_T_d_phi_d_T(indices) = ( V(:,:,ni)' * ...
      cds.dydt_ode(0,phases_T_i(:,i),parameters{:}) ) * delta_t(i) / period;
  end
 
  lhs_1_3 = V_T__b_g;
  
  for i=1:m % m == cds.nMeshIntervals
    indices           = (i-1)*basis_size + (1:basis_size);
    ni                = next_index_in_cycle(i,m);
    lhs_1_3(indices)  = lhs_1_3(indices) + V(:,:,ni)' * G_delta_q_gamma(:,i);
  end

  d_s_d_T      = 0;
  d_s_d_gamma  = 0;

  lhs_2_1 = [cds.previous_dydt_0' * V(:,:,1)    zeros(1,(m-1)*basis_size)];
  
  
  
  
  reduced_jacobian = [
    F_0_pp     V_T_d_phi_d_T     lhs_1_3     ;
    lhs_2_1    d_s_d_T           d_s_d_gamma + cds.previous_dydt_0' * delta_q_gamma(:,1);
  ];

end

  
function V = compute_subspace(i, period, parameters)
  global cds
  
  [eigenvectors, eigenvalues, no_convergence] = eigs( ...
    @(x) NewtonPicard.MultipleShooting.monodromy_map( ...
       i, x, period, parameters), ...
    cds.nphases, ...
    min(cds.nphases,cds.p + 1));


  if no_convergence
    V = [];
    fprintf(['Newton_Picard_Correction.m:', ...
      ' eigenvalues of monodromy matrix did not converge.\n'])
    return
  end

  eigenvalues = diag(eigenvalues);
  basis = zeros(cds.nphases, cds.p + 1);

  i = 1;

  while i <= cds.p
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
