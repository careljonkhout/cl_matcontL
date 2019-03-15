function v = find_tangent_vector(x)

  global cds

  active_par_val               = x(end);
  period                       = x(end-1);
  phases_0                     = x(1:end-2);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = active_par_val;
  parameters                   = num2cell(parameters);
  
  integration_opt = odeset(...
    'AbsTol',      1e-13,    ...
    'RelTol',      1e-13,    ...
    'BDF',         'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1,      ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  

  cds.cycle_trajectory = ode15s(...
    @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
    linspace(0, period, cds.nDiscretizationPoints), ...
    phases_0, integration_opt);

  phi = deval(cds.cycle_trajectory,period);

  
  if true || ~ isfield(cds, 'V')
    cds.V = compute_subspace(period, parameters);
    V = cds.V;
    basis_size = size(V,2);
  else
    cds.V = continue_subspace_with_convergence_criterium(period, parameters);
    V = cds.V;
    basis_size = size(V,2);
  end


  MV = zeros(cds.nphases,basis_size);
  
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
  dydt_monodromy_map = @(t, y) ...
    cds.jacobian_ode(t, deval(cds.cycle_trajectory,t), parameters{:}) * y;
  % The function monodromy_map cannot be used here, since it depends on
  % the global variable cds, and global variables are not copied so the
  % the workspace of the workers that parfor uses.
  for i=1:basis_size
    [~,trajectory] = ode15s(dydt_monodromy_map, [0 period], V(:,i), int_opt);
    MV(:,i) = trajectory(end,:)';
  end

  S = V'*MV;
  
  
  % the r in q_r means residual
  delta_q_r   = (eye(cds.nphases) - V*V')*(phi - phases_0);

  
  d_phi_d_gamma_val = d_phi_d_gamma(phases_0,period,parameters);
  
    delta_q_gamma = (eye(cds.nphases) - V*V') * d_phi_d_gamma_val;
  M_delta_q_gamma = monodromy_map(delta_q_gamma, period, parameters);
  
  
  d_phi_d_T = cds.dydt_ode(0,phi,parameters{:});
  

  d_s_d_T      = 0;
  d_s_d_gamma  = 0;

  
  left_hand_side = [
    S - eye(basis_size)         V'* d_phi_d_T     V' * (d_phi_d_gamma_val + M_delta_q_gamma);
    cds.previous_dydt_0' * V    d_s_d_T           d_s_d_gamma                               ;
  ];

  delta_p__delta_T_and_delta_gamma = null(left_hand_side);
  delta_p     = delta_p__delta_T_and_delta_gamma(1:end-2);
  delta_T     = delta_p__delta_T_and_delta_gamma(end-1);
  delta_gamma = delta_p__delta_T_and_delta_gamma(end);

  
  delta_q        = delta_q_r + delta_gamma * delta_q_gamma;
  delta_phases   = V * delta_p + delta_q;
  v              = [delta_phases; delta_T; delta_gamma];



function Mx = monodromy_map(phases_0, period, parameters)
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
    cds.jacobian_ode(t, deval(cds.cycle_trajectory,t), parameters{:}) * y;
  [~,trajectory] = ode15s(dydt_mon, [0 period], phases_0, int_opt);
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
 
function dphidp = d_phi_d_gamma(x, period, parameters)
  global cds
  ap = cds.ActiveParams;
  h = 1e-6;
  parameters{ap} = parameters{ap} - h;
  phi_1 = shoot(x, period, parameters);
  parameters{ap} = parameters{ap} + 2*h;
  phi_2 = shoot(x, period, parameters);
  dphidp = (phi_2 - phi_1)/h/2; 

  
function V = compute_subspace(period, parameters)
  global cds
  
  p = min([4 cds.nphases]);
  cds.p_extra = 2;
  cds.p = p;

  [eigenvectors, eigenvalues, no_convergence] = eigs( ...
    @(x) monodromy_map(x, period, parameters), ...
    cds.nphases, p + cds.p_extra);


  if no_convergence
    V = [];
    fprintf(['Newton_Picard_Correction.m:', ...
      ' eigenvalues of monodromy matrix did not converge.\n'])
    return
  end

  eigenvalues = diag(eigenvalues);
  basis = zeros(cds.nphases, p + cds.p_extra);
  cds.eigenvalues = eigenvalues;
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
