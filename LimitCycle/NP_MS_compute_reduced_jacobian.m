function [V, reduced_jacobian, delta_q_gamma, delta_q_r, G_delta_q_r, ...
          phases_0, phases_T_i, period, active_par_val] = ...
            NP_MS_compute_reduced_jacobian(x)

  global cds
  
  [phases_0, period, parameters] = ...
          NP_MS_extract_phases_period_and_parameters(x);
        
  cds.phases_0   = phases_0;
  active_par_val = parameters{cds.ActiveParams};
  m              = cds.n_mesh_intervals;
  mesh           = cds.mesh;
  
  % compute length of each mesh interval:
  delta_t               = period * diff(mesh);
   
  % when not using cvode, compute_cycle_parts must happen before computing
  % subspaces
  phases_T_i = NP_MS_compute_cycle_parts(x);

  V          = cell(m,1);
  V{1}       = compute_subspace(1, period, parameters);
  MV         = cell(m,1);
  M          = @NP_MS_monodromy_map; 
  
  for i = 1 : m % m == cds.n_mesh_intervals
    for j = 1 : size(V{i}, 2)
      MV{i}(:, j) = M(i, V{i}(:, j), delta_t(i), parameters);
    end
    if i == m
      break
    end
    
    V{i+1} = orth(MV{i});
    
    
    if size(V{i+1}, 2) == size(V{i}, 2); continue; end
    
    print_diag(1,'subspace is rank deficient. computing subspace using eigs\n');
    V{i+1} = compute_subspace(i + 1, period, parameters);
    
    if size(V{i+1},2) == size(V{i}, 2); continue; end
    
    cds.new_n_mesh_intervals = 1.5 * cds.n_mesh_intervals;
    % todo increase number of mesh intervals here
    % todo call "adapt" to remesh after increase of number of mesh intervals
     error('npms:wrong_basis_size', ...
         'subspace rank is still wrong. Abborting current continuation step\n');
   
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
  for i = 1 : m % note that: m == cds.n_mesh_intervals
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
  

  rhs_delta_q_gamma = zeros(cds.n_phases, m);
 
  V_T__b_g = zeros(reduced_jac_size, 1);
  
  row_offset = 0;
  for i = 1 : m
    b_g_i                  = NP_compute_d_phi_d_p(phases_0(:,i), ...
                                           delta_t(i), parameters);
    ni                     = next_index_in_cycle(i,m);
    indices                = row_offset + (1 : size(V{ni},2));
    V_T__b_g(indices)      = V{ni}' * b_g_i;
    
    rhs_delta_q_gamma(:,i) = b_g_i - V{ni} * V_T__b_g(indices);
    
    row_offset             = row_offset + size(V{ni},2);
  end
  
  
  % the r in q_r means residual
  
  rhs_delta_q_r = zeros(cds.n_phases,m);
  
  
  for i=1:m % m == cds.n_mesh_intervals
    ni                 = next_index_in_cycle(i,m);
    r                  = phases_T_i(:,i) - phases_0(:,ni);
    rhs_delta_q_r(:,i) = r - V{ni} * (V{ni}' * r);
  end
  
  solve_Q_system = @NP_MS_solve_Q_system;
  
  [delta_q_gamma, G_delta_q_gamma] = solve_Q_system(V, rhs_delta_q_gamma, ...
                                            delta_t, parameters);

  [delta_q_r    , G_delta_q_r    ] = solve_Q_system(V, rhs_delta_q_r, ...
                                            delta_t, parameters);
 
  jac_1_3 = V_T__b_g;
  
  row_offset = 0;
  for i = 1 : m % m == cds.n_mesh_intervals
    ni                = next_index_in_cycle(i,m);
    indices           = row_offset + (1 : size(V{ni}, 2));
    jac_1_3(indices)  = jac_1_3(indices) + V{ni}' * G_delta_q_gamma(:,i);
    row_offset        = row_offset + size(V{ni}, 2);
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

  
function V = compute_subspace_with_eigs(i, period, parameters)
  global cds
  
  n_eigenvalues = min(cds.n_phases, cds.preferred_basis_size + 1);

  [eigenvectors, eigenvalues] = eigs( ...
    @(x) NP_MS_full_monodromy_map(...
          i, x, period, parameters), ...
    cds.n_phases, ...
    n_eigenvalues);

  eigenvalues = diag(eigenvalues);
  basis = zeros(cds.n_phases, cds.preferred_basis_size + 1);

  j = 1;

  while j <= cds.preferred_basis_size
    basis(:, j) = real(eigenvectors(:,j));
    j = j + 1;
    if max(abs(imag(eigenvalues(j-1)))) > 1e-14
      basis(:, j) = imag(eigenvectors(:,j-1));
      j = j + 1;
    end
  end
  V = orth(basis(:,1:j-1));
end

function V = compute_subspace(i, period, parameters)
  global cds
  
  if ~ cds.using_cvode
    V = compute_subspace_with_eigs(i, period, parameters);
    return;
  end
  
  n_eigenvalues = min(cds.n_phases, cds.preferred_basis_size + 1);
  options.k = n_eigenvalues;
  options.blsz = 1;
  [Q, T] = ahbschur( ...
    @(X) multiple_monodromy_map_cvode(i, X, period, parameters), ...
    cds.n_phases, ...
    options);

  basis = zeros(cds.n_phases, cds.preferred_basis_size + 1);

  j = 1;

  while j <= cds.preferred_basis_size
    basis(:, j) = Q(:, j);
    j = j + 1;
    if T(j, j-1) ~= 0
      basis(:, j) = Q(:, j + 1);
      j = j + 1;
    end
  end
  V = basis(:,1:j-1);
end

function MX = multiple_monodromy_map(i, X, period, parameters)
  global cds contopts
  
  args.n_mesh_intervals = cds.n_mesh_intervals;
  args.mesh             = cds.mesh;
  args.phases_0         = cds.phases_0;
  args.integrator       = cds.integrator;
  args.i                = i;
  args.period           = period;
  args.parameters       = cell2mat(parameters);
  args.abs_tol          = contopts.integration_abs_tol;
  args.rel_tol          = contopts.integration_rel_tol;
  
  
  MX = zeros(size(X));
  for j = 1 : size(X,2)
    MX(:,j) = full_monodromy_parfor_enabled(X(:,j), args);
  end
end

function Mx = full_monodromy_parfor_enabled(x, args)
  Mx = x;
  delta_t = diff(args.mesh) * args.period;
  mesh_interval_sequence = [ args.i : args.n_mesh_intervals,  1 : args.i-1 ];
  for j = mesh_interval_sequence
    [~,~,Mx] = feval(args.integrator, ...
      't_values',                [0 delta_t(j)], ...
      'initial_point',           args.phases_0(:,j), ...
      'sensitivity_vector',      Mx, ...
      'ode_parameters',          args.parameters, ...
      'abs_tol',                 args.abs_tol, ...
      'rel_tol',                 args.rel_tol);
  end
end

function MX = multiple_monodromy_map_cvode(i, X, period, parameters)
  global cds contopts
  MX = X;
  delta_t = diff(cds.mesh) * period;
  mesh_interval_sequence = [ i : cds.n_mesh_intervals, 1 : i - 1];
  for j = mesh_interval_sequence
    [~, ~, MX] = feval(cds.integrator, ...
      't_values', [0 delta_t(j)], ...
      'initial_point', cds.phases_0(:, j), ...
      'sensitivity_vector', MX, ...
      'ode_parameters', cell2mat(parameters), ...
      'abs_tol', contopts.integration_abs_tol, ...
      'rel_tol', contopts.integration_rel_tol);
  end
  
end
                  
    



function index = next_index_in_cycle(i,m)
  index = i + 1;
  if index == m + 1
    index = 1;
  end
end
