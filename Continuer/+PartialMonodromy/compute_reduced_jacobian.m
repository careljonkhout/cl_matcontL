function [V, reduced_jacobian, ...
          phases_0, phi, period, active_par_val] = compute_reduced_jacobian(x)

  global cds contopts
  
  [phases_0,period,parameters] = ...
    NewtonPicard.SingleShooting.extract_phases_period_and_parameters(x);
  active_par_val = parameters{cds.ActiveParams};
  
  % Since the result of the next time-integration will be reused many times,
  % we set the tolerances a bit tighter than the rest.
  integration_opt = odeset(...
    'AbsTol',      contopts.integration_abs_tol / 10,    ...
    'RelTol',      contopts.integration_rel_tol / 10,    ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  

  cds.cycle_orbit = cds.integrator(...
    @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
    linspace(0, period, cds.nDiscretizationPoints), ...
    phases_0, integration_opt);

  phi = deval(cds.cycle_orbit,period);

  
  if ~ isfield(cds, 'V') || true
    V = orth(rand(length(phases_0),length(phases_0)/2));
    basis_size = size(V,2);
  end
  

  MV = zeros(cds.nphases,basis_size);
  
  int_opt = odeset(...
    'AbsTol',       contopts.integration_abs_tol,    ...
    'RelTol',       contopts.integration_rel_tol,    ...
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
  
  
  
  % the r in q_r means residual

  
  d_phi_d_T = cds.dydt_ode(0,phi,parameters{:});
  
  d_phi_d_gamma_val = d_phi_d_gamma(phases_0,period,parameters);
  d_s_d_T      = 0;
  d_s_d_gamma  = 0;

  
  reduced_jacobian = [
    V * MV'                     d_phi_d_T   d_phi_d_gamma_val;
    cds.previous_dydt_0'        d_s_d_T     d_s_d_gamma;
  ];

 

    
function x_end = shoot(x, period, parameters)
  global cds contopts
  f =@(t, y) cds.dydt_ode(t, y, parameters{:});
  integration_opt = odeset(...
    'AbsTol',      contopts.integration_abs_tol,    ...
    'RelTol',      contopts.integration_rel_tol,    ...
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

  
