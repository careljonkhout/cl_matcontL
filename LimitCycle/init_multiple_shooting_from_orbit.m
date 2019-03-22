%
% Initializes the global struct variable cds, and computes
% initial_continuation_data from an orbit that has converged to a limit cycle
% supplied by the user.
% 
% inputs:
% 
% orbit -- An orbit suplied by the user that has converged to a limit cycle.
% Only the last point of the orbit is used. The second dimension of the orbit
% array should correspond the dimensions of the system of ode's.
%
% nMeshPoints -- The number of points used in multiple shooting 
%
% odefile -- The odefile with the definition of the system of ODE's in standard
% matcontL format.
% 
% ode_parameters -- a cell array of the values of the parameters of the system
% of ODE's for which the orbit that converges to the limitcycle was found
%
% active_parameter_index -- the index ( in ode_parameters ) of the parameter
% that is varied in the continuation.
%
% time_integration_method -- one of the methods of the matlab ode suite, or a
% fully compatible alternative. (ode15s seems to works well for spatially
% discretized pde's)
%
% lower_bound_period -- a lower bound for the period of the cycle.
%
% upper_bound_period -- an upperbound for the period of the cycle.
%
% time_integration_options -- the options structure that are supplied to the
% time_integation_method. Use the matlab command odeset to generate this
% structure. See the documentation of odeset.
%
% poincare_tolerance -- The cycle is determined to be "closed", if the
% max-norm-distance from the starting point is less than poincare_tolerance (the
% end point of orbit). If one is trying to start a cycle that loops around twice
% before closing, one might want to set this to a small value, and supply very
% small tolerances to time_integration_options. Otherwise, don't set it too
% small (e.g. not smaller than 1e-2 for cycle with amplitudes between 1 and 10),
% or the closing of the cycle might not be detected.
function initial_continuation_data = init_multiple_shooting_from_orbit(...
    orbit, ...
    nMeshPoints, ...
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
  tangent_to_limitcycle    = dydt_ode(0,point_on_limitcycle, ode_parameters{:});
  time_integration_options = odeset(time_integration_options, ...
    'Events',       @returnToPlane, ...
    'Jacobian',     @(t,y) jacobian_ode(0, y, ode_parameters{:}));
  
  new_orbit = feval(time_integration_method, ...
    @(t,y) dydt_ode(0, y, ode_parameters{:}), ...
    linspace(0,upper_bound_period,cds.nDiscretizationPoints), ...
    point_on_limitcycle, ...
    time_integration_options); 
  
  
  period = new_orbit.x(end);
  initial_continuation_data = zeros(cds.nphases * nMeshPoints + 2,1);
  for i=0:nMeshPoints-1
    indices = (1:cds.nphases) + i * cds.nphases;
    initial_continuation_data(indices) = ...
      deval(new_orbit, i / nMeshPoints * period);
  end
  
  initial_continuation_data(end-1) = period;
  initial_continuation_data(end  ) = ode_parameters{active_parameter_index};
  
  cds.previous_phases = point_on_limitcycle;
  cds.previous_dydt_0 = tangent_to_limitcycle;
  cds.nMeshPoints = nMeshPoints;
  cds.dydt_ode        = dydt_ode;
  cds.jacobian_ode    = jacobian_ode;
  cds.integrator      = time_integration_method;
 
  function [value, isterminal, direction] = returnToPlane(t, x)
    % x and should be a column vector
    value = tangent_to_limitcycle'*(x-point_on_limitcycle);
    isterminal = t > lower_bound_period ...
      && max(abs(x-point_on_limitcycle)) < poincare_tolerance;
    direction = 1;
  end
end