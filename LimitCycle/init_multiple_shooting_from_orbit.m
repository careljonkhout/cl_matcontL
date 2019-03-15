function initial_continuation_data = init_multiple_shooting_from_orbit(...
    orbit, ...
    nShootingPoints, ...
    odefile, ...
    ode_parameters, ...
    active_parameter_index, ...
    time_integration_method, ...
    lower_bound_period, ...
    upper_bound_period, ...
    time_integration_options, ...
    poincare_tolerance)
  global cds
  
  if nargin <= 8
    time_integration_options = odeset(...
      'AbsTol',      1e-10,    ...
      'RelTol',      1e-10);
  end
  if nargin <= 9
    poincare_tolerance = 1e-2;
  end
  
  handles      = feval(odefile);
  dydt_ode     = handles{2};
  jacobian_ode = handles{3};

  point_on_limitcycle      = orbit(end,:)';
  tangent_to_limitcycle    = dydt_ode(0, point_on_limitcycle, ode_parameters{:});
  time_integration_options = odeset(time_integration_options, ...
    'Events',       @returnToPlane, ...
    'Jacobian',     @(t,y) jacobian_ode(0, y, ode_parameters{:}));
  
  new_orbit = feval(time_integration_method, ...
    @(t,y) dydt_ode(0, y, ode_parameters{:}), ...
    linspace(0,upper_bound_period,cds.nDiscretizationPoints), ...
    point_on_limitcycle, ...
    time_integration_options); 
  
  
  period = new_orbit.x(end);
  initial_continuation_data = zeros(cds.nphases * nShootingPoints + 2,1);
  for i=0:nShootingPoints-1
    indices = (1:cds.nphases) + i * cds.nphases;
    initial_continuation_data(indices) = ...
      deval(new_orbit, i / nShootingPoints * period);
  end
  
  initial_continuation_data(end-1) = period;
  initial_continuation_data(end  ) = ode_parameters{active_parameter_index};
  
  cds.previous_phases = point_on_limitcycle;
  cds.previous_dydt_0 = tangent_to_limitcycle;
  cds.nShootingPoints = nShootingPoints;
  cds.dydt_ode        = dydt_ode;
  cds.jacobian_ode    = jacobian_ode;
  
  function [value, isterminal, direction] = returnToPlane(t, x)
    % x and should be a column vector
    value = tangent_to_limitcycle'*(x-point_on_limitcycle);
    isterminal = t > lower_bound_period ...
      && max(abs(x-point_on_limitcycle)) < poincare_tolerance;
    direction = 1;
  end
end