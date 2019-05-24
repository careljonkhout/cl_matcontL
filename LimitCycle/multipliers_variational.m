function multipliers = multipliers_variational(x)
  global lds cds;
  print_diag(3,'computing multipliers\n');
  points_on_cycle              = x(1:end-2);
  points_on_cycle              = reshape(points_on_cycle, ...
                                         cds.nphases, length(lds.finemsh));
  period                       = x(end);
  active_parameter_value       = x(end-1);
  parameters                   = num2cell(lds.P0);
  parameters{lds.ActiveParams} = active_parameter_value;
  
  
 [~, multiplier_matrix, no_convergence] = eigs( ...
      @(x) monodromy_map(x, points_on_cycle, period, parameters), ...
      cds.nphases, cds.nphases);
  multipliers = diag(multiplier_matrix);

  if no_convergence
    print_diag(2, 'multipliers did not converge\n');
    return
  end
  
  distance_to_one = abs(multipliers - 1);
  accuracy        = min(distance_to_one);
  print_diag(0,'deviation of trivial multiplier: %.2e\n', accuracy);
  print_diag(2, 'multiplier: %.12f\n', multipliers);
end

function Mx  = monodromy_map(phases_0, ...
  points_on_cycle, time_interval, parameters)
  global cds lds contopts
  
  jacobian = @(t,y) cds.jacobian_ode(t, ...
                interp1(lds.finemsh, points_on_cycle', t, 'spline'), ...
                parameters{:});

  int_opt = odeset(...
    'AbsTol',       contopts.integration_abs_tol,    ...
    'RelTol',       contopts.integration_rel_tol,    ...
    'Jacobian',     jacobian ...
  );

  dydt_mon = @(t, y) jacobian(t, y) * y;
  
  [~,orbit] = cds.integrator(dydt_mon, [0 time_interval], phases_0, int_opt);
  
  Mx = orbit(end,:)';
end