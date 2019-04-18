function multipliers = compute_multipliers(x, nMults_to_compute)
  global cds contopts;
  print_diag(3,'computing multipliers\n');
  [~, period, parameters] = ...
        NewtonPicard.MultipleShooting.extract_phases_period_and_parameters(x);
  
  NewtonPicard.MultipleShooting.compute_stiched_orbit(x, ...
                                             contopts.multipliers_rel_tol, ...
                                             contopts.multipliers_abs_tol);
  
  [~, multiplier_matrix, no_convergence] = eigs( ...
     @(x) monodromy_map(x, period, parameters), cds.nphases, nMults_to_compute);
  multipliers = diag(multiplier_matrix);
  
  if no_convergence
    print_diag(2, 'multipliers did not converge\n');
    return
  end
  
  print_diag(2, multipliers2str(multipliers));
end

function Mx  = monodromy_map(phases_0, delta_t, parameters)
  global cds contopts
  int_opt = odeset(...
    'AbsTol',       contopts.multipliers_abs_tol,    ...
    'RelTol',       contopts.multipliers_rel_tol,    ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode, ...
                      t, interp1(cds.t_cycle,cds.y_cycle,t,'spline'), ...
                      parameters{:}) ...
  );

  dydt_mon = @(t, y) ...
    cds.jacobian_ode(t,interp1(cds.t_cycle,cds.y_cycle,t,'spline'), ...
      parameters{:}) * y;
  
  [~,orbit] = cds.integrator(dydt_mon, [0 delta_t], phases_0, int_opt);
  
  Mx = orbit(end,:)';
end