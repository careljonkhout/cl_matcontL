function Mx  = NP_MS_full_monodromy_map(i, x, period, parameters)
  global cds contopts
  if ~ contopts.monodromy_by_finite_differences
    Mx = x;
    delta_t = diff(cds.mesh) * period;
    mesh_interval_sequence = [ i : cds.n_mesh_intervals,  1 : i-1 ];
    for j = mesh_interval_sequence
      Mx = NP_MS_monodromy_map(j, Mx, delta_t(j), ...
                                  parameters, ...
                                  contopts.multipliers_abs_tol, ...
                                  contopts.multipliers_rel_tol);
    end
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
    ff = @(t, x1_and_x2) [f(0, x1_and_x2(1:cds.n_phases    ))
                          f(0, x1_and_x2(  cds.n_phases+1:end))];
    [~, orbit] = ode15s(ff, ...
      [0 period], ...
      [x_cycle - h * x; x_cycle+h * x], ...
      integration_opt);
    
    phi_x1__and__phi_x2 = orbit(end,:)';
    phi_x1 = phi_x1__and__phi_x2(1:cds.n_phases);
    phi_x2 = phi_x1__and__phi_x2(cds.n_phases+1:end);
    
    Mx = (phi_x2 - phi_x1)/h/2; 
  end
end