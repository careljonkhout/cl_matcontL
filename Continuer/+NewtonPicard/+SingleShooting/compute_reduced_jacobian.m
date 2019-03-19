% Computes the Jacobian matrix with subspace reduced applied to the phase space
% also returns additional values that were used in the computation of the
% reduces Jacobian matrix, which are needed in do_one_correction.m and/or 
% find_tangent_vector.m. The subspace used is the subspace associated to the
% largest eigenvalues (in modulus) of the mondromy matrix.

function [V, reduced_jacobian, delta_q_gamma, delta_q_r, M_delta_q_r, ...
          phases_0, phi, period, active_par_val] = compute_reduced_jacobian(x)

  global cds
  
  [phases_0,period,parameters] = ...
    NewtonPicard.SingleShooting.extract_phases_period_and_parameters(x);
  active_par_val = parameters{cds.ActiveParams};
  
  integration_opt = odeset(...
    'AbsTol',      1e-10,    ...
    'RelTol',      1e-10,    ...
    'BDF',         'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1,      ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  

  cds.cycle_trajectory = cds.integrator(...
    @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
    linspace(0, period, cds.nDiscretizationPoints), ...
    phases_0, integration_opt);

  phi = deval(cds.cycle_trajectory,period);

  
  if ~ isfield(cds, 'V')
    cds.V = NewtonPicard.SingleShooting.compute_subspace(period, parameters);
    V = cds.V;
    basis_size = size(V,2);
  else
    cds.V = NewtonPicard.SingleShooting.continue_subspace(period, parameters);
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
    [~,trajectory] = cds.integrator(dydt_monodromy_map, [0 period], V(:,i), int_opt);
    MV(:,i) = trajectory(end,:)';
  end

  S = V'*MV;
  
  
  % the r in q_r means residual
  [delta_q_r,     M_delta_q_r] = ...
    NewtonPicard.SingleShooting.solve_Q_system(V, phi - phases_0, ...
    period, parameters);

  d_phi_d_gamma_val = d_phi_d_gamma(phases_0,period,parameters);
  [delta_q_gamma, M_delta_q_gamma] = ...
    NewtonPicard.SingleShooting.solve_Q_system(V, d_phi_d_gamma_val, ...
    period, parameters);
  

  
  d_phi_d_T = cds.dydt_ode(0,phi,parameters{:});
  

  d_s_d_T      = 0;
  d_s_d_gamma  = 0;

  
  reduced_jacobian = [
    S - eye(basis_size)         V'* d_phi_d_T     V' * (d_phi_d_gamma_val + M_delta_q_gamma);
    cds.previous_dydt_0' * V    d_s_d_T           d_s_d_gamma                               ;
  ];

 

    
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
  [~, trajectory] = cds.integrator(f, [0 period], x, integration_opt);
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

  
