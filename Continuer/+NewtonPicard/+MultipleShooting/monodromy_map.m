%
% i:                index of shooting-point from where to start
% phases_0:         vector to which the monodromy map is applied
% time_interval:    length of the time interval for time integration
% parameters:       cell array of parameters for the jacobian of the ode
function Mx  = monodromy_map(i, phases_0, time_interval, parameters)
  global cds contopts
  int_opt = odeset(...
    'AbsTol',       contopts.monodromy_map_abs_tol,    ...
    'RelTol',       contopts.monodromy_map_rel_tol,    ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode, ...
                      t, deval(cds.orbits(i),t), parameters{:}) ...
  );

  dydt_mon = @(t, y) ...
    cds.jacobian_ode(t, deval(cds.orbits(i), t), parameters{:}) * y;
  
  [~,orbit] = cds.integrator(...
    dydt_mon, [0 time_interval], phases_0, int_opt);
  
  Mx = orbit(end,:)';
end