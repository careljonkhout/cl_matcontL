function multipliers = compute_multipliers(x)
  global cds;
  print_diag(3,'computing multipliers\n');
  [phases, period, parameters] = ...
    NewtonPicard.SingleShooting.extract_phases_period_and_parameters(x);
  
  integration_opt = odeset(...
    'AbsTol',      1e-13,    ...
    'RelTol',      1e-13,    ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  
  cds.cycle_trajectory = ode15s(...
    @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
    linspace(0, period, cds.nDiscretizationPoints), ...
    phases, integration_opt);
  
  [~, multiplier_matrix, no_convergence] = eigs( ...
    @(x) NewtonPicard.SingleShooting.monodromy_map(x, period, parameters), ...
    cds.nphases, cds.p);
  multipliers = diag(multiplier_matrix);
  
  if no_convergence
    print_diag(2, 'multipliers did not converge\n');
  else
    print_diag(2, 'multiplier: %.12f\n', multipliers);
  end
end