% based on algorithm outlined in paragraph 6.2 of (bibtex citation follows)
% @phdthesis{lust-phd,
%	  author={Lust, Kurt},
%	  title={Numerical bifurcation analysis 
%        of periodic solutions of partial differential equations},
%	  school={K.U.Leuven},
%	  year={1997},
% }
% most variable names are derived from variable names in \cite{lust-phd}.
function x = Newton_Picard_Multiple_Shooting(x0,x,v0)

  global cds
  [phases_0,period,parameters] = getComponents(x);
  active_par_val               = x(end);
  
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
   
  if ~ isfield(cds, 'V') || true
    for i=2:m % m == cds.nShootingPoints
      V(:,:,i) = compute_subspace(i, period, parameters);
    end
  else
    %for i=1:m
    %  cds.V = continue_subspace_with_convergence_criterium(i,period,parameters);
    %end
    %V = cds.V;
    %basis_size = size(V,2);
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
  
  
  [delta_q_gamma,G_delta_q_gamma] = solve_q_systems(V, ...
    rhs_delta_q_gamma, partial_period, parameters, cds.q_systems_tolerance);

  [delta_q_r    ,G_delta_q_r]     = solve_q_systems(V, ...
    rhs_delta_q_r, partial_period, parameters, cds.q_systems_tolerance);
  
  
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
  lhs_3_1 = zeros(1,basis_size * m);
  for i=1:m % m == cds.nShootingPoints
    indices_lhs_3_1 = (i-1)*basis_size  + (1:basis_size);
    indices_v0      = (i-1)*cds.nphases + (1:cds.nphases);
    lhs_3_1(indices_lhs_3_1) = v0(indices_v0)' * V(:,:,i);
  end
  
  left_hand_side = [
    F_0_pp     V_T_d_phi_d_T     lhs_1_3     ;
    lhs_2_1    d_s_d_T           d_s_d_gamma ;
    lhs_3_1    v0(end-1)         v0(end)     ;
  ];


  rhs_1 = zeros(basis_size*m,1);
  for i=1:m
    indices        = (i-1) * basis_size + (1:basis_size);
    ni             = next_index_in_cycle(i,m);
    rhs_1(indices) = V(:,:,ni)'*(phi(:,i) - phases_0(:,ni) + G_delta_q_r(:,i));
  end
 
  
  right_hand_side = - [
    rhs_1;
    
    cds.previous_dydt_0' * ...
      (phases_0(:,1) - cds.previous_phases  + delta_q_r(:,1));
      
    v0'*(x-x0)  + v0(1:end-2)' * reshape(delta_q_r,numel(delta_q_r),1);
  ];

  delta_p__delta_T_and_delta_gamma = left_hand_side \ right_hand_side;
  delta_p     = delta_p__delta_T_and_delta_gamma(1:end-2);
  delta_T     = delta_p__delta_T_and_delta_gamma(end-1);
  delta_gamma = delta_p__delta_T_and_delta_gamma(end);

  V_delta_p = zeros(cds.nphases*m,1);
  for i=1:m
    indices1 = (i-1) * cds.nphases + (1:cds.nphases);
    indices2 = (i-1) * basis_size  + (1:basis_size );
    V_delta_p(indices1) = V(:,:,i) * delta_p(indices2);
  end
  phases_0 = reshape(phases_0,numel(phases_0),1);
  delta_q         = delta_q_r + delta_gamma * delta_q_gamma;
  phases_0        = phases_0 + V_delta_p + reshape(delta_q,numel(delta_q),1);
  period          = period + delta_T;
  active_par_val  = active_par_val + delta_gamma;
  x = [phases_0; period; active_par_val];
  %v = find_tangent_vector(phases_0, period, parameters, V);
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


function Mx = monodromy_map(i, phases_0, time_interval, parameters)
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

  while i <= cds.p + 1 - 2
    basis(:, i) = real(eigenvectors(:,i));
    i = i + 1;
    if abs(imag(eigenvalues(i-1))) > 0
      basis(:, i) = imag(eigenvectors(:,i-1));
      i = i + 1;
    end
  end
  V = orth(basis(:,1:i));
end

% based on page 283 of (bibtex citation follows)
% @incollection{lust2000,
%	title={Computation and bifurcation analysis of
%                          periodic solutions of large-scale systems},
%	author={Lust, Kurt and Roose, Dirk},
% booktitle={Numerical methods for bifurcation problems 
%                          and large-scale dynamical systems},
%	pages={265--301},
%	year={2000},
%	publisher={Springer}
% }
function V = continue_subspace_with_convergence_criterium(period, parameters) %#ok<DEFNU>
  global cds contopts
  V_extended = [cds.V rand(cds.nphases, cds.p + cds.p_extra - size(cds.V,2))];
  %V_extended = [rand(cds.nphases, cds.p + cds.p_extra)];
  p = cds.p;
  int_opt = odeset(...
    'AbsTol',       1e-10,    ...
    'RelTol',       1e-10,    ...
    'BDF',          'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1, ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode, ...
                      t, deval(cds.cycle_trajectory,t), parameters{:}) ...
  );                
  dydt_partial_monodromy_map = @(t, y) ...
    cds.jacobian_ode(t, deval(cds.cycle_trajectory,t), parameters{:}) * y;

  
  p_eff = 0;
  W = V_extended;
  iteration = 0;
  while p_eff < p && iteration < contopts.NewtonPicardMaxSubspaceIterations
    iteration = iteration + 1;
    if  iteration > 3
      fprintf('subspace iteration %d\n',iteration);
      V_extended(:,p_eff+1:end) = orth(W(:,p_eff+1:end));
      %V_extended = orth(W);
    end
    
    try
      % one can try parfor here, but errors might occur (see catch block)
      for i=p_eff+1:size(V_extended,2)
        % The function partial_monodromy_map cannot be used here,
        % since it depends on the global variable cds, and global variables 
        % are not copied to the workspaces of the workers that parfor uses.
        [~,trajectory] = ode15s(...
          dydt_partial_monodromy_map, [0 period], V_extended(:,i), int_opt);
        W(:,i) = trajectory(end,:)';
      end
    catch error
      if strcmp(error.identifier,'MATLAB:remoteparfor:AllParforWorkersAborted')
        % Something went wrong with the parfor workers.
        % We try again with ordinary for.
        fprintf('Parfor aborted, retrying with ordinary for.\n');
        for i=p_eff+1:size(V_extended,2)
          [~,trajectory] = ode15s(...
            dydt_partial_monodromy_map, [0 period], V_extended(:,i), int_opt);
          W(:,i) = trajectory(end,:)';
        end
      else
        % in case of some other error, we want to know about it
        rethrow(error)
      end
    end
        
      
    U = V_extended'*W;
    [Y,S] = schur(U);
    % order schur vectors
    [~,I] = sort(abs(ordeig(S)));
    I(I) = 1:length(I);
    [Y,S] = ordschur(Y,S,I);
    V_extended = V_extended * Y;
    W = W * Y;
    k = size(V_extended,2) - 1;
    % note: svds(A,1) is a built-in matlab function
    % that computes the largest singular value of A    
    while k >= 1 ...
        && (S(k+1,k) ~= 0 || S(k+1,k) == 0 ...
        && svds(W(:,1:k) - V_extended(:,1:k) * S(1:k,1:k), 1) ...
                 >= contopts.NewtonPicardBasisTolerance )
      k = k - 1;

    end
    p_eff = k;
    %fprintf('p_eff: %d subspace norm at k=p_eff: %.10f\n', p_eff, ...
    %  svds(W(:,1:p_eff) - V_extended(:,1:p_eff) * S(1:p_eff,1:p_eff), 1));
  end
  V = V_extended(:,1:p_eff);
end

  
% inputs (Lust)
% based on algorithm 6.2 on pages 197 and 198 of (bibtex citation follows)
% @phdthesis{lust-phd,
%	  author={Lust, Kurt},
%	  title={Numerical bifurcation analysis
%          of periodic solutions of partial differential equations},
%	  school={K.U.Leuven},
%	  year={1997},
% }
% inputs:
% - starting values delta_q_i, i=0..m-1, each in \R^N (N is number of spatial 
%   dimensions of ode
% - optionally G_i(delta_q_i) i=0..m-1 each in \R^n
% - convergence thresholds eps_i 
% - bases V for projectors
% - right hand sides rhs
% - routine "partial_monodromy_map" to compute G delta_q_i
%
% P will be the projectors onto the small subspaces
%
function [delta_q, G_delta_q] = ...
             solve_q_systems(V, rhs, partial_period, parameters, tolerance)
  global cds;
  m = cds.nShootingPoints;
  for i=1:m
    ni = next_index_in_cycle(i,m);
    rhs(:,i) = rhs(:,i) - V(:,:,ni) * V(:,:,ni)' * rhs(:,i);
  end
  
%   if max(max(abs(rhs))) < eps
%     delta_q   = zeros(cds.nphases,m);
%     G_delta_q = zeros(cds.nphases,m);
%     return
%   end
  G_delta_q  = zeros(cds.nphases,m);
  delta_q    = zeros(cds.nphases,m);
  for iteration_number = 1:15
   
    for i=2:m
      delta_q(:,i) = G_delta_q(:,i-1) + rhs(:,i-1);
      delta_q(:,i) = delta_q(:,i) - V(:,:,i) * V(:,:,i)' * delta_q(:,i);
      G_delta_q(:,i) = monodromy_map(i, ...
        delta_q(:,i), partial_period, parameters);
    end
    condensed_residual = G_delta_q(:,m) + rhs(:,m);
    condensed_residual = condensed_residual ...
                         - V(:,:,1) * V(:,:,1)' * condensed_residual;
    if max(max(abs(delta_q(:,1) - condensed_residual))) < tolerance
      break
    else
      delta_q(:,1) = condensed_residual;
      G_delta_q(:,1) = monodromy_map(1, delta_q(:,1), ...
        partial_period, parameters);
    end
  end
  print_diag(5,'did %d iterations in solve_q_systems\n',iteration_number);
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
  
    
