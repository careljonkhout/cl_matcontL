% Computes the Jacobian matrix with subspace reduced applied to the phase space
% also returns additional values that were used in the computation of the
% reduces Jacobian matrix, which are needed in do_one_correction.m and/or 
% find_tangent_vector.m. The subspace used is the subspace associated to the
% largest eigenvalues (in modulus) of the mondromy matrix.

function [V, reduced_jacobian, delta_q_gamma, delta_q_r, M_delta_q_r, ...
          phases_0, phi, period, active_par_val] = compute_reduced_jacobian(x)

  global cds contopts
  
  [phases_0,period,parameters] = ...
    NewtonPicard.SingleShooting.extract_phases_period_and_parameters(x);
  active_par_val = parameters{cds.ActiveParams};
  
  integration_opt = odeset(...
    'AbsTol',      contopts.orbit_abs_tol,    ...
    'RelTol',      contopts.orbit_rel_tol,    ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  

  cds.cycle_orbit = cds.integrator(...
    @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
    linspace(0, period, cds.nDiscretizationPoints), ...
    phases_0, integration_opt);

  phi = deval(cds.cycle_orbit,period);

  
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
    'AbsTol',       contopts.MV_abs_tol,    ...
    'RelTol',       contopts.MV_rel_tol,    ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode, ...
                      t, deval(cds.cycle_orbit,t), parameters{:}) ...
  );                
  dydt_monodromy_map = @(t, y) ...
    cds.jacobian_ode(t, deval(cds.cycle_orbit,t), parameters{:}) * y;
  % The function monodromy_map cannot be used here, since it depends on
  % the global variable cds, and global variables are not copied so the
  % the workspace of the workers that parfor uses.
  integrator = cds.integrator;
  if contopts.contL_ParallelComputing
    try 
      parfor i=1:basis_size
        % The function monodromy_map cannot be used here, since it depends on
        % the global variable cds, and global variables are not copied so the
        % the workspace of the workers that parfor uses.
        [~,orbit] = feval(integrator, ...
          dydt_monodromy_map, [0 period], V(:,i), int_opt);

        MV(:,i) = orbit(end,:)';
      end
      parfor_failed = false;
    catch error
      if (strcmp(error.identifier,'MATLAB:remoteparfor:AllParforWorkersAborted'))
        % Something went wrong with the parfor workers.
        % We try again with ordinary for.
        fprintf('Parfor aborted, retrying with ordinary for.\n');
        parfor_failed = true;
      else
        % in case of some other error, we want to know about it
        rethrow(error)
      end
    end
  end
  
  if (~ contopts.contL_ParallelComputing) || parfor_failed
    for i=1:basis_size
      [~,orbit] = integrator(dydt_monodromy_map, [0 period], V(:,i), int_opt);
      MV(:,i) = orbit(end,:)';
    end
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
  global cds contopts
  f =@(t, y) cds.dydt_ode(t, y, parameters{:});
  integration_opt = odeset(...
    'AbsTol',      contopts.shoot_abs_tol,    ...
    'RelTol',      contopts.shoot_rel_tol,    ...
    'BDF',         'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1,      ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  [~, orbit] = cds.integrator(f, [0 period], x, integration_opt);
  x_end = orbit(end,:)';
 
function dphidp = d_phi_d_gamma(x, period, parameters)
  global cds
  ap = cds.ActiveParams;
  h = 1e-6;
  parameters{ap} = parameters{ap} - h;
  phi_1 = shoot(x, period, parameters);
  parameters{ap} = parameters{ap} + 2*h;
  phi_2 = shoot(x, period, parameters);
  dphidp = (phi_2 - phi_1)/h/2; 

  
