function multipliers = compute_multipliers(x, nMults_to_compute)
  global cds contopts;
  print_diag(3,'computing multipliers\n');
  [phases, period, parameters] = ...
    NewtonPicard.SingleShooting.extract_phases_period_and_parameters(x);
  
  integration_opt = odeset(...
    'AbsTol',      contopts.integration_abs_tol,    ...
    'RelTol',      contopts.integration_rel_tol,    ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  
  cds.cycle_orbit = cds.integrator(...
    @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
    linspace(0, period, cds.nDiscretizationPoints), ...
    phases(:,1), integration_opt);
  
  [~, multiplier_matrix, no_convergence] = eigs( ...
    @(x) NewtonPicard.SingleShooting.monodromy_map(x, period, parameters), ...
    cds.nphases, nMults_to_compute);
  multipliers = diag(multiplier_matrix);
  
  if no_convergence
    print_diag(2, 'multipliers did not converge\n');
    return
  end
  print_diag(4, 'nMults_to_compute: %d\n', nMults_to_compute);
  print_diag(4, 'multipliers computed: %d\n', length(multipliers));
  distance_to_one = abs(multipliers - 1);
  accuracy        = min(distance_to_one);
  if nargin == 1
    print_diag(0,'deviation of trivial multiplier: %.2e\n', accuracy);
  end
  print_diag(2, ['multipliers:\n' multipliers2str(multipliers)]);
  
end