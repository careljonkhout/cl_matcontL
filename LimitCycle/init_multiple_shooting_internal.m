% Intended to be called by other functions. To initialize a cycle continuation
% using multiple shooting call one of these two functions:
%
% init_multiple_shooting_find_stable_cycle
% init_multiple_shooting_from_hopf

function initial_continuation_data = init_multiple_shooting_internal(in)

  global cds
  cds = [];
  
  handles      = feval(in.odefile);
  dydt_ode     = handles{2};
  jacobian_ode = handles{3};
  cds.n_phases = length(in.point_on_limitcycle);
  norm_of_gap  = Inf;
  
  tangent_to_limitcycle =  ...
          dydt_ode(0,in.point_on_limitcycle, in.ode_parameters{:});
  in.time_integration_options = odeset(in.time_integration_options, ...
          'Events',       @returnToPlane);
  
  if ~ isempty(jacobian_ode)
    in.time_integration_options = odeset(in.time_integration_options, ...
      'Jacobian',     @(t,y) jacobian_ode(0, y, in.ode_parameters{:}));
  end
  
  using_cvode = endsWith(func2str(in.time_integration_method), 'cvode');

  
  if using_cvode
    [orbit_t,orbit_y] = feval(in.time_integration_method, ...
      't_values',                linspace(0, in.upper_bound_period, 100), ...
      'initial_point',           in.point_on_limitcycle, ...
      'cycle_detection',         true, ...
      'lower_bound_period',      in.lower_bound_period, ...
      'point_on_limitcycle',     in.point_on_limitcycle, ...
      'ode_parameters',          cell2mat(in.ode_parameters), ...
      'abs_tol',                 in.time_integration_options.AbsTol, ...
      'rel_tol',                 in.time_integration_options.RelTol);

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
  end
  
  period   = orbit_t(end);
  
  matlab_gui_running = usejava('jvm');
  
  if in.show_plots && matlab_gui_running
    my_figure = figure;
    plot_t = linspace(0, period, in.n_interpolated_points);
    transformed_orbit = in.plot_transformation(orbit_y) ...
                      - in.plot_transformation(orbit_y(1,:));
    plot_y = interp1(orbit_t, transformed_orbit, plot_t, in.interpolation);
    plot(plot_t', plot_y);
    xlabel('t')
    ylabel('deviation form initial value')
    disp(['Now showing plot from  t = time_to_converge_to_cycle  to  ' ...
                                 't = time_to_converge_to_cycle + period']);
    my_pause();
    if isvalid(my_figure)
      close(my_figure.Number)
    end
  end
  
  if ~ using_cvode
    fprintf('max-norm of gap in cycle: %.4e\n', norm_of_gap)
  end    
  
  
  period   = orbit_t(end);

  
  initial_continuation_data = zeros(cds.n_phases * in.n_mesh_intervals + 2, 1);
  if using_cvode
    m = in.n_mesh_intervals;
    t_values = linspace(0, period * (m - 1) / m, m);
    [~, orbit_y] = feval(in.time_integration_method, ...
      't_values',                t_values, ...
      'initial_point',           in.point_on_limitcycle, ...
      'ode_parameters',          cell2mat(in.ode_parameters), ...
      'abs_tol',                 in.time_integration_options.AbsTol, ...
      'rel_tol',                 in.time_integration_options.RelTol);
    
    orbit_y = orbit_y';
    initial_continuation_data(1:end-2) = orbit_y(:);
  else
    orbit = feval(in.time_integration_method, ...
        @(t,y) dydt_ode(0, y, in.ode_parameters{:}), ...
        [0 period], ...
        in.point_on_limitcycle, ...
        in.time_integration_options);
      
    
    for i = 0 : in.n_mesh_intervals - 1
      indices = (1:cds.n_phases) + i * cds.n_phases;
      initial_continuation_data(indices) = ...
        deval(orbit, i / in.n_mesh_intervals * period);
    end
  end
  
  initial_continuation_data(end-1) = period;
  initial_continuation_data(end  ) = ...
                                   in.ode_parameters{in.active_parameter_index};
  
  cds.n_mesh_intervals     = in.n_mesh_intervals;
  cds.new_n_mesh_intervals = in.n_mesh_intervals;
  cds.mesh                 = linspace(0, 1, in.n_mesh_intervals + 1);
  cds.probfile             = in.odefile;
  cds.options.PartitionMonodromy = cds.n_phases > 20;
  cds.nap                  = 1;
  cds.ndim                 = cds.n_phases * cds.n_mesh_intervals + 2;
  cds.usernorm             = [];
  cds.ncoo                 = cds.n_phases * cds.n_mesh_intervals + 1;
  cds.ActiveParams         = in.active_parameter_index;
  cds.P0                   = cell2mat(in.ode_parameters);
  cds.previous_phases      = in.point_on_limitcycle(:);
  cds.previous_dydt_0      = tangent_to_limitcycle(:);
  cds.dydt_ode             = dydt_ode;
  cds.jacobian_ode         = jacobian_ode;
  cds.jacobian_p_ode       = handles{4};
  cds.integrator           = in.time_integration_method;
  cds.preferred_basis_size  = in.subspace_size;
  cds.p                    = in.subspace_size;
  cds.mv_count             = 0;
  cds.curve                = @multiple_shooting;
  cds.using_cvode          = using_cvode;
 
  function [value, isterminal, direction] = returnToPlane(t, x)
    % x and should be a column vector
    value = tangent_to_limitcycle'*(x-in.point_on_limitcycle);
    norm_of_gap = max(abs(x-in.point_on_limitcycle));
    isterminal = t > in.lower_bound_period && ...
                     norm_of_gap < in.poincare_tolerance;
    direction = 1;
  end
end