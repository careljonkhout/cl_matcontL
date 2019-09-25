function multipliers = compute_multipliers(x, nMults_to_compute)
  global cds contopts;
  print_diag(3,'computing multipliers\n');
  [~, period, parameters] = ...
        NewtonPicard.MultipleShooting.extract_phases_period_and_parameters(x);
  
  if ( ~ contopts.monodromy_by_finite_differences) && ( ~ cds.using_cvode )
    NewtonPicard.MultipleShooting.compute_stiched_orbit(x, ...
                                             contopts.multipliers_rel_tol, ...
                                             contopts.multipliers_abs_tol);
  end
  
  if cds.using_cvode
    M = @(x) ...
          NewtonPicard.MultipleShooting.monodromy_map(1, x, period, parameters);
  else
    M = @(x) monodromy_map(x, period, parameters);
  end
  
  nMults_to_compute = min(nMults_to_compute, cds.nphases);
  
  [~, multiplier_matrix, no_convergence] = eigs(...
                                             M, cds.nphases, nMults_to_compute);
  multipliers = diag(multiplier_matrix);
  
  if no_convergence
    print_diag(2, 'multipliers did not converge\n');
    return
  end
  
  print_diag(1, multipliers2str(multipliers));
end

function Mx  = monodromy_map(phases_0, delta_t, parameters)
  global cds contopts
  if ~ contopts.monodromy_by_finite_differences
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
  else
    % Below an alternative method of computing the action of the monodromy
    % matrix is implemented. Here, the action of the monodromy matrix is
    % computed by finite differences. Seems to be faster than the method above.
    % The error of the trivial multiplier is much larger when using finite
    % differences.
 
    
    integration_opt = odeset(...
      'AbsTol',       1e-13, ...
      'RelTol',       1e-13  ...
    ); 
    h = 5e-5;
    x_cycle = deval(cds.orbits(1),0);
    f  = @(t, x) cds.dydt_ode(0,x,parameters{:});
    ff = @(t, x1_and_x2) [f(0, x1_and_x2(1:cds.nphases    ))
                          f(0, x1_and_x2(  cds.nphases+1:end))];
    [~, orbit] = ode15s(ff, ...
      [0 delta_t], ...
      [x_cycle - h * phases_0; x_cycle+h * phases_0], ...
      integration_opt);
    
    phi_x1__and__phi_x2 = orbit(end,:)';
    phi_x1 = phi_x1__and__phi_x2(1:cds.nphases);
    phi_x2 = phi_x1__and__phi_x2(cds.nphases+1:end);
    
    Mx = (phi_x2 - phi_x1)/h/2; 
  end
end