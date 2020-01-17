function phases_T_i = NP_MS_compute_cycle_parts(x)
  global cds contopts
  m = cds.n_mesh_intervals;

  [phases_0, period, parameters] = ...
          NP_MS_extract_phases_period_and_parameters(x);
        
  delta_t = period * diff(cds.mesh);
  
  if ~ cds.using_cvode
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
  end
  
  phases_T_i = zeros(cds.n_phases, m);
  
  if cds.using_cvode
    common_args.abs_tol    = contopts.integration_abs_tol;
    common_args.rel_tol    = contopts.integration_rel_tol;
    common_args.integrator = cds.integrator;
    common_args.parameters = cell2mat(parameters);
    if contopts.contL_ParallelComputing
      % note: although parfor is used here, it does not lead to significant
      % speedups, since this loop is not the most time consuming loop in
      % multiple shooting continuations.
      parfor i = 1 : m
        phases_T_i(:,i) = shoot_cvode(phases_0(:,i), delta_t(i), common_args);
      end
    else
      for i = 1 : m
        phases_T_i(:,i) = shoot_cvode(phases_0(:,i), delta_t(i), common_args);
      end
    end
  else
    for i = 1 : m
      cds.orbits(i) = cds.integrator(...
        @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
        [0 delta_t(i) * 1.1], ...
        phases_0(:,i), integration_opt);

      phases_T_i(:,i) = deval(cds.orbits(i), delta_t(i));
    end
  end
end

function phases_T_i = shoot_cvode(phases_0, delta_t, args)
  [~, y] = feval(args.integrator, ...
      'initial_point',   phases_0, ...
      't_values',        [0 delta_t], ...
      'ode_parameters',  args.parameters, ...
      'abs_tol',         args.abs_tol, ...
      'rel_tol',         args.rel_tol);
    phases_T_i = y(end,:)';
end
    