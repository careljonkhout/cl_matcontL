function point = Newton_Picard_Corrections(x)
  
  global cds contopts
 
  
  curve_function_norm = max(abs(feval(cds.curve_func,x)));
  corrections = 0;
  while (~ (curve_function_norm < contopts.FunTolerance)) ...
      && corrections < contopts.MaxCorrIters
    
    corrections = corrections + 1;
    period = x(end-1);
    fprintf('function_norm: %.8e period: %.8e corrections: %d\n', ...
      curve_function_norm, period, corrections);
    x = Newton_Picard_Correction(x);
    if isempty(x)
      return;
    end

    if period < 0
      break
    end
    curve_function_norm = norm(feval(cds.curve_func,x));
  end
  if curve_function_norm < contopts.FunTolerance
    point.R = curve_function_norm; 
    point.x = x;
    point.iters = corrections;
  else
    point = [];
    return
  end
  
  

function [x,v] = Newton_Picard_Correction(x0)

  global cds

  active_par_val               = x0(end);
  period                       = x0(end-1);
  phases_0                     = x0(1:end-2);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = active_par_val;
  parameters                   = num2cell(parameters);
  
  integration_opt = odeset(...
    'AbsTol',      1e-10,    ...
    'RelTol',      1e-10,    ...
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

  p = 10;

  [eigenvectors, eigenvalues, no_convergence] = eigs( ...
    @(x) monodromy_map(x, period, parameters), ...
    cds.nphases, p);
%'StartVector' — Initial starting vector 

  if no_convergence
    x = [];
    fprintf(['Newton_Picard_Correction.m:', ...
      ' eigenvalues of monodromy matrix did not converge.\n'])
    return
  end

  eigenvalues = diag(eigenvalues);
  basis = zeros(cds.nphases,p-2);

  basis_index = 0;
  i = 0;

  while basis_index <= p - 2
    i = i + 1;
    basis_index = basis_index + 1;
    basis(:, basis_index) = real(eigenvectors(:,i));
    if abs(imag(eigenvalues(i))) > 0
      basis_index = basis_index + 1;
      i = i + 1;
      basis(:, basis_index) = imag(eigenvectors(:,i));
    end
  end

  basis_size = basis_index;
  V = orth(basis);

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
  parfor i=1:basis_size
    [~,trajectory] = ode15s(dydt_monodromy_map, [0 period], V(:,i), int_opt);
    MV(:,i) = trajectory(end,:)';
  end

  S = V'*MV;
  delta_q = (eye(cds.nphases) - V*V')*(phi - phases_0);
  M_delta_q = monodromy_map(delta_q, period, parameters);

  d_phi_d_T = cds.dydt_ode(0,phi,parameters{:});

  d_s_d_T = cds.dydt_0' * d_phi_d_T;

  lhs = [S - eye(basis_size) V'* d_phi_d_T];
  lhs = [lhs; cds.dydt_0' * V d_s_d_T];

  rhs = - V'*(phi - phases_0 + M_delta_q);
  rhs = [rhs; -cds.dydt_0' * (phi - phases_0 + delta_q)];

  xi_and_dT = lhs\rhs;
  xi = xi_and_dT(1:end-1);
  dT = xi_and_dT(end);

  phases = phases_0 + V*xi + delta_q;
  period = period + dT;
  x = [phases; period; active_par_val];
  %v = find_tangent_vector(phases_0, period, parameters, V);



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
    

  
% not used 
function v = find_tangent_vector(phases_0,period,parameters,V) %#ok<DEFNU>
  global cds
  approximate_monodromy = V*V';
  jacobian              = [approximate_monodromy-eye(cds.nphases); cds.dydt_0'];
  phases_end            = shoot(phases_0, period, parameters);
  d_phi_d_T             = cds.dydt_ode(0,phases_end,parameters{:});
  d_s_d_T               = cds.dydt_0' * d_phi_d_T;
  jacobian              = [jacobian [d_phi_d_T; d_s_d_T]];
  dphidp                = d_phi_d_p(phases_0, period, parameters);
  dsdp                  = cds.dydt_0'*dphidp;
  jacobian              = [jacobian [dphidp; dsdp]];
  v                     = null(jacobian);
  
% only used by find_tangent_vector
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
 
% only used by find_tangent vector
function dphidp = d_phi_d_p(x, period, parameters)
  global cds
  ap = cds.ActiveParams;
  h = 1e-6;
  parameters{ap} = parameters{ap} - h;
  phi_1 = shoot(x, period, parameters);
  parameters{ap} = parameters{ap} + 2*h;
  phi_2 = shoot(x, period, parameters);
  dphidp = (phi_2 - phi_1)/h/2; 

  
