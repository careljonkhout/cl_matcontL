function phases_T_i = compute_cycle_parts(x)
  global cds contopts
  m = cds.n_mesh_intervals;

  [phases_0, period, parameters] = ...
          NewtonPicard.MultipleShooting.extract_phases_period_and_parameters(x);
        
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
  for i=1:m
    print_diag(6,'orbit %d\n',i)
    if cds.using_cvode
      [~, y] = cds.integrator( ...
        'initial_point',   phases_0(:,i), ...
        't_values',        [0 delta_t(i)], ...
        'ode_parameters',  cell2mat(parameters), ...
        'abs_tol',         contopts.integration_abs_tol, ...
        'rel_tol',         contopts.integration_rel_tol);
      phases_T_i(:,i) = y(end,:)';
    else
      cds.orbits(i) = cds.integrator(...
        @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
        [0 delta_t(i) * 1.1], ...
        phases_0(:,i), integration_opt);

      phases_T_i(:,i) = deval(cds.orbits(i), delta_t(i));
    end
  end
  
end