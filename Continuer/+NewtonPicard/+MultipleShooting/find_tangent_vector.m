% based on algorithm outlined in paragraph 6.2 of (bibtex citation follows)
% @phdthesis{lust-phd,
%	  author={Lust, Kurt},
%	  title={Numerical bifurcation analysis 
%        of periodic solutions of partial differential equations},
%	  school={K.U.Leuven},
%	  year={1997},
% }
% most variable names are derived from variable names in \cite{lust-phd}.
function v = find_tangent_vector(x)

  global cds
  [phases_0,period,parameters] = getComponents(x);
  
  integration_opt = odeset(...
    'AbsTol',      1e-10,    ...
    'RelTol',      1e-10,    ...
    'BDF',         'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1,      ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  
  m = cds.nShootingPoints;
  phi = zeros(cds.nphases,m);

  for i=1:m
    cds.trajectories(i) = ode15s(...
      @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
      linspace(0, period, cds.nDiscretizationPoints), ...
      phases_0(:,i), integration_opt);
    
    phi(:,i) = deval(cds.trajectories(i), period/m);
  end
  
  cds.p = cds.preferred_basis_size;
  
  V1             = compute_subspace(1,period, parameters);
  basis_size     = size(V1,2);
  cds.basis_size = basis_size;
  V              = zeros(cds.nphases,basis_size,m);
  V(:,:,1)       = V1;
   
  if ~ isfield(cds, 'V')
    for i=2:m % m == cds.nShootingPoints

      %compute_subspace(i, period, parameters);


     for j = 1:size(V,2)
       V(:,j,i) = monodromy_map(i-1, V(:,j,i-1), period/m, parameters);
     end
     V(:,:,i) = orth(V(:,:,i));

    end
    cds.V = V;
  else
    V(:,:,1) = NewtonPicard.MultipleShooting.continue_subspace( ...
        1, phases_0(:,1), period, parameters);
    for i=2:m
      for j = 1:size(V,2)
        V(:,j,i) = monodromy_map(i-1, V(:,j,i-1), period/m, parameters);
      end
      V(:,:,i) = orth(V(:,:,i));
    end
    cds.V = V;
  end


  MV = zeros(cds.nphases,basis_size,cds.nShootingPoints);
 
  int_opt = odeset(...
    'AbsTol',       1e-10,    ...
    'RelTol',       1e-10,    ...
    'BDF',          'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1);

  % the variable name F_0_pp corresponds with L^0_pp in Lust
  % we compute the number of nonzero's (nnz):
  nnz = m* basis_size * basis_size; % blocks on diagonal;
  nnz = nnz + m * basis_size;       % identity matrices on superdiagonal
  F_0_pp = spalloc(m*basis_size,m*basis_size,nnz);
  
  for i=1:m % note that: m == cds.nShootingPoints
    int_opt = odeset(int_opt, ...
      'Jacobian', @(t,y) feval(cds.jacobian_ode, ...
                    t, deval(cds.trajectories(i),t), parameters{:}));
    dydt_partial_mon = @(t, y) ...
      cds.jacobian_ode(t,deval(cds.trajectories(i),t),parameters{:}) * y;
    % The function partial_monodromy_map cannot be used here, in case we want to 
    % use parfor, since it depends on
    % the global variable cds, and global variables are not copied so the
    % the workspace of the workers that parfor uses.
    for j=1:basis_size
      [~,trajectory] = ...
        ode15s(dydt_partial_mon, [0 period/m], V(:,j,i), int_opt);
      MV(:,j,i) = trajectory(end,:)';
    end

    indices = basis_size*(i-1) + (1:basis_size);
    ni = next_index_in_cycle(i,m);
    F_0_pp(indices,indices)= V(:,:,ni)' * MV(:,:,i); %#ok<SPRIX>
    
    if i < m
      F_0_pp(indices, indices + basis_size) = - eye(basis_size); %#ok<SPRIX>
    else
      F_0_pp(indices, 1:basis_size)         = - eye(basis_size); %#ok<SPRIX>
    end
  end
  
  
  
  V_T__b_T = zeros(basis_size * m,1); % V_p^T b_T      in \cite{lust-phd}
  
  indices = 1:basis_size;
  for i=1:m
    b_T_i             = cds.dydt_ode(0, phi(:,i), parameters{:}) / m;
    V_T__b_T(indices) = V(:,:,next_index_in_cycle(i,m))' * b_T_i;
    indices           = indices + basis_size;
  end
  
 
  rhs_delta_q_gamma = zeros(cds.nphases,m);
 
  V_T__b_g = zeros(basis_size * m,1);
 
  
  indices = 1:basis_size;
  for i=1:m
    b_g_i                  = d_phi_d_gamma(phases_0(:,i), period/m, parameters);
    
    
    ni                     = next_index_in_cycle(i,m);
    V_T__b_g(indices)      = V(:,:,ni)' * b_g_i;
    
    rhs_delta_q_gamma(:,i) = b_g_i - V(:,:,ni) * V_T__b_g(indices);
    indices                = indices + basis_size;
  end
  
  
  % the r in q_r means residual
  
  rhs_delta_q_r = zeros(cds.nphases,m);
  
  for i=1:m % m == cds.nShootingPoints
    ni                 = next_index_in_cycle(i,m);
    r                  = phi(:,i) - phases_0(:,ni);
    rhs_delta_q_r(:,i) = r - V(:,:,ni) * V(:,:,ni)' * r;
  end

  %q_systems_tolerance = 1e-6;
  
  partial_period = period / m; % m == cds.nShootingpoints
  
  
  [delta_q_gamma,G_delta_q_gamma] = ...
    NewtonPicard.MultipleShooting.solve_Q_system(V, ...
    rhs_delta_q_gamma, partial_period, parameters);

  [delta_q_r    ,~]     = ...
    NewtonPicard.MultipleShooting.solve_Q_system(V, ...
    rhs_delta_q_r, partial_period, parameters);
  
  
  V_T_d_phi_d_T = zeros(basis_size*m,1);
  for i=1:m % m == cds.nShootingPoints
    indices = (i-1)*basis_size + (1:basis_size);
    ni = next_index_in_cycle(i,m);
    V_T_d_phi_d_T(indices) ...
      = V(:,:,ni)' * cds.dydt_ode(0,phi(:,i),parameters{:}) /m;
  end
 
  lhs_1_3 = V_T__b_g;
  
  for i=1:m % m == cds.nShootingPoints
    indices           = (i-1)*basis_size + (1:basis_size);
    ni                = next_index_in_cycle(i,m);
    lhs_1_3(indices)  = lhs_1_3(indices) + V(:,:,ni)' * G_delta_q_gamma(:,i);
  end

  d_s_d_T      = 0;
  d_s_d_gamma  = 0;

  lhs_2_1 = [cds.previous_dydt_0' * V(:,:,1)    zeros(1,(m-1)*basis_size)];
  
  left_hand_side = [
    F_0_pp     V_T_d_phi_d_T     lhs_1_3     ;
    lhs_2_1    d_s_d_T           d_s_d_gamma ;
  ];

  delta_p__delta_T_and_delta_gamma = null(full(left_hand_side));
  delta_p     = delta_p__delta_T_and_delta_gamma(1:end-2);
  delta_T     = delta_p__delta_T_and_delta_gamma(end-1);
  delta_gamma = delta_p__delta_T_and_delta_gamma(end);

  V_delta_p = zeros(cds.nphases*m,1);
  for i=1:m
    indices1 = (i-1) * cds.nphases + (1:cds.nphases);
    indices2 = (i-1) * basis_size  + (1:basis_size );
    V_delta_p(indices1) = V(:,:,i) * delta_p(indices2);
  end
  
  delta_q         = delta_q_r + delta_gamma * delta_q_gamma;
  delta_phases        = V_delta_p + reshape(delta_q,numel(delta_q),1);
  v = [delta_phases; delta_T; delta_gamma];
end

% extracts 
% - y ( the current approximation of points on the cycle )
% - period
% - parameters ( of the ode system in which cycles are continued )
% from the continuation state vector x.
% The non-active parameters, i.e. the parameters that do not change during
% the continuation are extracted from the global struct cds
% (i.e.) curve description structure.
% The parameters are returned as a cell array, so that that can be passed to
% cds.dydt_ode in an syntactically elegant manner.

function [y,period,parameters] = getComponents(x)
  global cds
  y                            = x(1:cds.nphases*cds.nShootingPoints);
  y                            = reshape(y,cds.nphases,cds.nShootingPoints);
  period                       = x(end-1);
  parameter_value              = x(end);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = parameter_value;
  parameters                   = num2cell(parameters);
end

%
% i:                index of shooting-point from where to start
% phases_0:         vector to which the monodromy map is applied
% time_interval:    length of the time interval for time integration
% parameters:       cell array of parameters for the jacobian of the ode
function Mx  = monodromy_map(i, phases_0, time_interval, parameters)
  global cds
  int_opt = odeset(...
    'AbsTol',       1e-10,    ...
    'RelTol',       1e-10,    ...
    'BDF',          'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1, ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode, ...
                      t, deval(cds.trajectories(i),t), parameters{:}) ...
  );                
  dydt_mon = @(t, y) ...
    cds.jacobian_ode(t, deval(cds.trajectories(i), t), parameters{:}) * y;
  [~,trajectory] = ode15s(dydt_mon, [0 time_interval], phases_0, int_opt);
  Mx = trajectory(end,:)';
end
    
function x_end = shoot(x, period, parameters)
  global cds
  f =@(t, y) cds.dydt_ode(t, y, parameters{:});
  integration_opt = odeset(...
    'AbsTol',      1e-10,    ...
    'RelTol',      1e-10,    ...
    'BDF',         'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1,      ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  [~, trajectory] = ode15s(f, [0 period], x, integration_opt);
  x_end = trajectory(end,:)';
end

 
function dphidp = d_phi_d_gamma(x, partial_period, parameters)
  global cds
  ap = cds.ActiveParams;
  h = 1e-6;
  parameters{ap} = parameters{ap} - h;
  phi_1 = shoot(x, partial_period, parameters);
  parameters{ap} = parameters{ap} + 2*h;
  phi_2 = shoot(x, partial_period, parameters);
  dphidp = (phi_2 - phi_1)/h/2;
end

  
function V = compute_subspace(i, period, parameters)
  global cds
  
  [eigenvectors, eigenvalues, no_convergence] = eigs( ...
    @(x) monodromy_map(i, x, period, parameters), ...
    cds.nphases, ...
    cds.p + 1);


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

function index = previous_index_in_cycle(i,m) %#ok<DEFNU>
  index = i - 1;
  if index == 0
    index = m;
  end
end
  
    
