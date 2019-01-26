function point = Newton_Picard_Multiple_Shooting(x,v)
  
  global cds contopts
 
  curve_function_norm = max(abs(feval(cds.curve_func,x)));
  corrections = 0;
  while (~ (curve_function_norm < contopts.FunTolerance)) ...
      && corrections < contopts.MaxCorrIters

    period = x(end-1);
    fprintf('function_norm: %.8e period: %.8e corrections: %d\n', ...
      curve_function_norm, period, corrections);
    corrections = corrections + 1;
    x = Newton_Picard_Correction(x,v);
    if isempty(x)
      break;
    end

    if period < 0
      fprintf('period less than zero, correction aborted\n');
      break
    end
    new_curve_function_norm = max(abs(feval(cds.curve_func,x)));
    if (new_curve_function_norm > curve_function_norm) 
      fprintf('function_norm: %.8e period: %.8e corrections: %d\n', ...
      new_curve_function_norm, period, corrections);
      fprintf('Curve function norm is increasing. Aborting Corrections.\n');
      break;
    end
    curve_function_norm = new_curve_function_norm;
  end
  if curve_function_norm < contopts.FunTolerance
    point.R = curve_function_norm; 
    point.x = x;
    point.v = v;
    point.iters = corrections;
  else
    point = [];
    return
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
  
  

% based on algorithm outlined in paragraph 6.2 of (bibtex citation follows)
% @phdthesis{lust-phd,
%	  author={Lust, Kurt},
%	  title={Numerical bifurcation analysis 
%        of periodic solutions of partial differential equations},
%	  school={K.U.Leuven},
%	  year={1997},
% }
% most variable names are derived from variable names in \cite{lust-phd}.
function x = Newton_Picard_Correction(x0,v0)

  global cds
  [phases_0,period,parameters] = getComponents(x0);
  active_par_val               = x0(end);
  
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
    cds.partial_trajectories(i) = ode15s(...
      @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
      linspace(0, period/m, cds.nDiscretizationPoints), ...
      phases_0(:,i), integration_opt);
    
    phi(:,i) = deval(cds.partial_trajectories(i), period/m);
  end

  
  if ~ isfield(cds, 'V') && true
    for i=1:m
      V(:,:,i) = compute_subspace(i, period, parameters);
    end
    basis_size = size(V,2);
  else
    for i=1:m
      cds.V = continue_subspace_with_convergence_criterium(i,period,parameters);
    end
    V = cds.V;
    basis_size = size(V,2);
  end


  MV = zeros(cds.nphases,basis_size,cds.nShootingPoints);
  
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
  
  %S = zeros(basis_size,basis_size,m);
  % the variable name F_0_pp corresponds with L^0_pp in Lust 
  F_0_pp = spalloc(m*basis_size,m*basis_size,m*basis_size*basis_size);
  
  for i=1:m
    dydt_partial_mon = @(t, y) ...
      cds.jacobian_ode(t,deval(cds.partial_trajectory(i),t),parameters{:}) * y;
    % The function partial_monodromy_map cannot be used here, in case we want to 
    % use parfor, since it depends on
    % the global variable cds, and global variables are not copied so the
    % the workspace of the workers that parfor uses.
    for i=1:basis_size
      [~,trajectory] = ode15s(dydt_partial_mon, [0 period], V(:,i), int_opt);
      MV(:,i,i) = trajectory(end,:)';
    end
    % S(:,:,i) = V(:,:,i)'*MV(:,:,i);  
    indices = basis_size*(i-1);
    F_0_pp(indices,indices)= V(:,:,i)' * MV(:,:,i); %#ok<SPRIX>
    if i < m
      F_0_pp(indices, indices + basis_size) = - eye(basis_size); %#ok<SPRIX>
    else
      F_0_pp(indices, 1:basis_size)       = - eye(basis_size); %#ok<SPRIX>
    end
  end
  
  
  
  V_T__b_T = zeros(basis_size * m); % V_p^T b_T      in \cite{lust-phd}
  b_T_i = cds.dydt_ode(0, phi(m), parameters{:}) / m;
  V_T__b_T(1:basis_size) = V(:,:,1)' * b_T_i;
  
  for i=2:m
    b_T_i             = cds.dydt_ode(0, phi(i), parameters{:}) / m;
    indices           = (i-1)*basis_size + 1:basis_size;
    V_T__b_T(indices) = V(:,:,i) * b_T_i;
  end
  
  V_T__b_g = zeros(basis_size * m); % V_p^T b_\gamma in \cite{lust-phd}
  b_g_i = d_phi_d_gamma(phases_0(:,m), period/m, parameters);
  V_T__b_g(1:basis_size) = V(:,:,1)' * b_g_i;
  
  for i=2:m
    b_g_i             = d_phi_d_gamma(phases_0(:,i-1), period/m, parameters);
    indices           = (i-1)*basis_size + 1:basis_size;
    V_T__b_g(indices) = V(:,:,i) * b_g_i;
  end
  
  
  % the r in q_r means residual
  

  
  d_phi_d_gamma_val = d_phi_d_gamma(phases_0,period,parameters);
  
    delta_q_gamma = (eye(cds.nphases) - V*V') * d_phi_d_gamma_val;
  M_delta_q_gamma = partial_monodromy_map(delta_q_gamma, period, parameters);
  
  
  d_phi_d_T = cds.dydt_ode(0,phi,parameters{:});
  

  d_s_d_T      = 0;
  d_s_d_gamma  = 0;

  
  left_hand_side = [
    S - eye(basis_size)         V'* d_phi_d_T     V' * (d_phi_d_gamma_val + M_delta_q_gamma);
    cds.previous_dydt_0' * V    d_s_d_T           d_s_d_gamma                               ;
    v0(1:end-2)' * V            v0(end-1)         v0(end)                                   ;
  ];

  right_hand_side = [
    - V'                   * (phi - phases_0                + M_delta_q_r);
    - cds.previous_dydt_0' * (phases_0 - cds.previous_phases  + delta_q_r);
    - v0(1:end-2)'         * (phi - phases_0                  + delta_q_r);
  ];%                          

  delta_p__delta_T_and_delta_gamma = left_hand_side \ right_hand_side;
  delta_p     = delta_p__delta_T_and_delta_gamma(1:end-2);
  delta_T     = delta_p__delta_T_and_delta_gamma(end-1);
  delta_gamma = delta_p__delta_T_and_delta_gamma(end);

  
  delta_q         = delta_q_r + delta_gamma * delta_q_gamma;
  phases_0          = phases_0 + V * delta_p + delta_q;
  period          = period + delta_T;
  active_par_val  = active_par_val + delta_gamma;
  x = [phases_0; period; active_par_val];
  %v = find_tangent_vector(phases_0, period, parameters, V);



function Mx = partial_monodromy_map(i, phases_0, partial_period, parameters)
  global cds
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
  dydt_mon = @(t, y) ...
    cds.jacobian_ode(t, deval(cds.partial_trajectory(i)), parameters{:}) * y;
  [~,trajectory] = ode15s(dydt_mon, [0 partial_period], phases_0, int_opt);
  Mx = trajectory(end,:)';
    
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
  phi_1 = shoot(x, period, parameters);
  parameters{ap} = parameters{ap} + 2*h;
  phi_2 = shoot(x, period, parameters);
  dphidp = (phi_2 - phi_1)/h/2; 

  
function V = compute_subspace(i, period, parameters)
  global cds
  
  p = min([10 cds.nphases]);
  cds.p_extra = 2;
  cds.p = p;

  [eigenvectors, eigenvalues, no_convergence] = eigs( ...
    @(x) partial_monodromy_map(i, x, period, parameters), ...
    cds.nphases, p + cds.p_extra);


  if no_convergence
    V = [];
    fprintf(['Newton_Picard_Correction.m:', ...
      ' eigenvalues of monodromy matrix did not converge.\n'])
    return
  end

  eigenvalues = diag(eigenvalues);
  basis = zeros(cds.nphases, p + cds.p_extra);

  i = 0;

  while i <= p + cds.p_extra - 1
    i = i + 1;
    basis(:, i) = real(eigenvectors(:,i));
    if abs(imag(eigenvalues(i))) > 0
      i = i + 1;
      basis(:, i) = imag(eigenvectors(:,i-1));
    end
  end
  V = orth(basis(:,1:i));

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
function V = continue_subspace_with_convergence_criterium(period, parameters)
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
                    solve_q_systems(delta_q_0, P, rhs, period, parameters, eps)
  global cds;
  N = size(V_p,1);
  m = cds.nShootingPoints;
  G_delta_q_0 = zeros(size(delta_q_0));
  G_delta_q_0(:,m) = ...
    partial_monodromy_map(1, delta_q_0(:,m  ), period/m, parameters);
  for i=2:m
    G_delta_q_0(:,i) = ...
      partial_monodromy_map(i, delta_q_0(:,i-1), period/m, parameters);
  end
  for i=1:m
    rhs(:,i) = rhs(:,i) + G_delta_q_0(:,i) + delta_q_0(:,i-1);
    rhs(:,i) = rhs(:,i) - P(:,:,i) * rhs(:,i);
  end
  if max(max(abs(rhs))) < eps
    delta_q = delta_q_0;
    G_delta_q = G_delta_q_0;
    return
  end
  i = 0;
  G_delta_q = zeros(size(delta_q_0));
  delta_q   = zeros(size(delta_q_0));
  done = false;
  while ~ done
    % according to Lust there is a potential optimization here
    % that is not implemented in this code for now
    G_delta_q(:,1) = partial_monodromy_map(...
                      1, delta_q(:,m), period/m, parameters);
    for i=1:m-1
      delta_q(:,i) = G_delta_q(:,i) + rhs(:,i);
      delta_q(:,i) = delta_q(:,i) - P * delta_q(:,i);
      G_delta_q(:,i+1) = partial_monodromy_map(...
                          i+1, delta_q(:,i), period/m, parameters);
    end
    v = rhs(:,m) + G_delta_q(:,m) - P * G_delta_q(:,m) - delta_q(:,m);
    done = max(max(abs(v))) < eps;
  end
  delta_q   =   delta_q +   delta_q_0;
  G_delta_q = G_delta_q + G_delta_q_0;
  
  

