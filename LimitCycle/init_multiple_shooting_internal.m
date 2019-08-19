% Intended to be called by other functions. To initialize a cycle continuation
% using multiple shooting call one of these two functions:
%
% init_multiple_shooting_find_stable_cycle
% init_multiple_shooting_from_hopf

function initial_continuation_data = init_multiple_shooting_internal(in)
  clear global
  global cds
  
  handles      = feval(in.odefile);
  dydt_ode     = handles{2};
  jacobian_ode = handles{3};
  cds.nphases  = length(in.point_on_limitcycle);
  
  tangent_to_limitcycle =  ...
    dydt_ode(0,in.point_on_limitcycle, in.ode_parameters{:});
  in.time_integration_options = odeset(in.time_integration_options, ...
    'Events',       @returnToPlane);
  
  if ~ isempty(jacobian_ode)
    in.time_integration_options = odeset(in.time_integration_options, ...
      'Jacobian',     @(t,y) jacobian_ode(0, y, in.ode_parameters{:}));
  end
  
  new_orbit = feval(in.time_integration_method, ...
    @(t,y) dydt_ode(0, y, in.ode_parameters{:}), ...
    [0 in.upper_bound_period], ...
    in.point_on_limitcycle, ...
    in.time_integration_options); 
  
  
  period                    = new_orbit.x(end);
  
  if in.show_plots
    trail_solution_t = linspace(0, period, 10000);
    trail_solution_x = deval(new_orbit, trail_solution_t);
    my_figure = figure;
    plot(trail_solution_t, trail_solution_x-new_orbit.y(:,1))
    xlabel('t')
    ylabel('deviation form initial value')
    disp(['Now showing plot from t=time_to_converge_to_cycle to ' ...
                                       't=time_to_converge_to_cycle + period']);
    disp('Press a key to continue')
    pause
    if isvalid(my_figure)
      close(my_figure.Number)
    end
  end
  
  initial_continuation_data = zeros(cds.nphases * in.nMeshIntervals + 2, 1);
  for i=0:in.nMeshIntervals-1
    indices = (1:cds.nphases) + i * cds.nphases;
    initial_continuation_data(indices) = ...
      deval(new_orbit, i / in.nMeshIntervals * period);
  end
  
  initial_continuation_data(end-1) = period;
  initial_continuation_data(end  ) = ...
    in.ode_parameters{in.active_parameter_index};
  
  cds.using_cvode     = false; % todo: add cvode support
  cds.nMeshIntervals  = in.nMeshIntervals;
  cds.mesh            = linspace(0, 1, in.nMeshIntervals + 1);
  cds.probfile        = in.odefile;
  cds.options.PartitionMonodromy = cds.nphases > 20;
  cds.nap             = 1;
  cds.ndim            = cds.nphases * cds.nMeshIntervals + 2;
  cds.usernorm        = [];
  cds.ncoo            = cds.nphases * cds.nMeshIntervals + 1;
  cds.ActiveParams    = in.active_parameter_index;
  cds.P0              = cell2mat(in.ode_parameters);
  cds.previous_phases = in.point_on_limitcycle;
  cds.previous_dydt_0 = tangent_to_limitcycle;
  cds.dydt_ode        = dydt_ode;
  cds.jacobian_ode    = jacobian_ode;
  cds.integrator      = in.time_integration_method;
  cds.preferred_basis_size  = in.subspace_size;
  cds.p               = in.subspace_size;
  cds.mv_count        = 0;
  cds.curve           = @multiple_shooting;
 
  function [value, isterminal, direction] = returnToPlane(t, x)
    % x and should be a column vector
    value = tangent_to_limitcycle'*(x-in.point_on_limitcycle);
    isterminal = t > in.lower_bound_period ...
      && max(abs(x-in.point_on_limitcycle)) < in.poincare_tolerance;
    direction = 1;
  end
end