% Intended to be called by other functions. To initialize a cycle continuation
% using single shooting call one of these two functions:
%
% init_single_shooting_find_stable_cycle
% init_single_shooting_from_hopf

function initial_continuation_data = init_single_shooting_internal(in)
  clear global  
  global cds
  
  
  handles      = feval(in.odefile);
  dydt_ode     = handles{2};
  jacobian_ode = handles{3};
  cds.nphases  = length(in.point_on_limitcycle);
  norm_of_gap  = Inf;
  
 

  tangent_to_limitcycle ...   
    = dydt_ode(0,in.point_on_limitcycle, in.ode_parameters{:});
  
  using_cvode = endsWith(func2str(in.time_integration_method), 'cvode');
  
  if using_cvode
    [orbit_t,orbit_y] = feval(in.time_integration_method, ...
      't_values',                linspace(0,in.upper_bound_period,2000), ...
      'initial_point',           in.point_on_limitcycle, ...
      'cycle_detection',         true, ...
      'lower_bound_period',      in.lower_bound_period, ...
      'point_on_limitcycle',     in.point_on_limitcycle, ...
      'ode_parameters',          cell2mat(in.ode_parameters), ...
      'abs_tol',                 in.time_integration_options.AbsTol, ...
      'rel_tol',                 in.time_integration_options.RelTol);

    if in.show_plots
      my_figure = figure;
      plot(orbit_t, orbit_y - orbit_y(1,:));
    end
  else
    in.time_integration_options = odeset(in.time_integration_options, ...
      'Events',       @returnToPlane);
    if ~ isempty(jacobian_ode)
      in.time_integration_options = odeset(in.time_integration_options, ...
        'Jacobian',     @(t,y) jacobian_ode(0, y, in.ode_parameters{:}));
    end

    [orbit_t, orbit_y] = feval(in.time_integration_method, ...
      @(t,y) dydt_ode(0, y, in.ode_parameters{:}), ...
      [0 in.upper_bound_period], ...
      in.point_on_limitcycle, ...
      in.time_integration_options);
    
    my_figure = figure;
    plot(orbit_t, orbit_y - orbit_y(1,:))
  end
  
  if in.show_plots
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
  
  if ~ using_cvode
    fprintf('max-norm of gap in cycle: %.4e\n', norm_of_gap)
  end
  
  
  period                    = orbit_t(end);
  initial_continuation_data = zeros(cds.nphases + 2, 1);
  initial_continuation_data(1:cds.nphases) = in.point_on_limitcycle;
  
  
  initial_continuation_data(end-1) = period;
  initial_continuation_data(end) = in.ode_parameters{in.active_parameter_index};
  
  cds.probfile        = in.odefile;
  cds.options.PartitionMonodromy = cds.nphases > 20;
  cds.nap             = 1;
  cds.ndim            = cds.nphases + 2;
  cds.usernorm        = [];
  cds.ncoo            = cds.nphases + 1;
  cds.ActiveParams    = in.active_parameter_index;
  cds.P0              = cell2mat(in.ode_parameters);
  cds.previous_phases = in.point_on_limitcycle(:);
  cds.previous_dydt_0 = tangent_to_limitcycle(:);
  cds.dydt_ode        = dydt_ode;
  cds.jacobian_ode    = jacobian_ode;
  cds.jacobian_p_ode  = handles{4};
  cds.integrator      = in.time_integration_method;
  cds.preferred_basis_size  = in.subspace_size;
  cds.p               = in.subspace_size;
  cds.mv_count        = 0;
  cds.curve           = @single_shooting;
  cds.using_cvode     = using_cvode;
 
  function [value, isterminal, direction] = returnToPlane(t, x)
    % x and should be a column vector
    value = tangent_to_limitcycle'*(x-in.point_on_limitcycle);
    norm_of_gap = max(abs(x-in.point_on_limitcycle));
    isterminal = t > in.lower_bound_period && ...
                     norm_of_gap < in.poincare_tolerance;
    direction = 1;
  end
end