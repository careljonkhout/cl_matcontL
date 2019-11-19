function point_on_cycle = converge_to_cycle(in)
    
  global cds;
  
  handles       = feval(in.odefile);
  dydt_ode      = handles{2};
  jacobian_ode  = handles{3};
  cds.n_phases  = length(in.initial_point);
  using_cvode   = endsWith(func2str(in.time_integration_method), 'cvode') || ...
                  endsWith(func2str(in.time_integration_method), 'arkode') ;
  
  if ~ isempty(jacobian_ode) && ~ using_cvode
    in.time_integration_options = odeset(in.time_integration_options, ...
                   'Jacobian', @(t,y) jacobian_ode(0, y, in.ode_parameters{:}));
  end
  
  orbit_to_cycle_t = linspace(0, in.time_to_converge_to_cycle, ...
                                 in.n_computed_points);
  
  if using_cvode
    [~, orbit_to_cycle_y]= feval(in.time_integration_method, ...
      't_values',                orbit_to_cycle_t, ...
      'initial_point',           in.initial_point, ...
      'ode_parameters',          cell2mat(in.ode_parameters), ...
      'abs_tol',                 in.time_integration_options.AbsTol, ...
      'rel_tol',                 in.time_integration_options.RelTol, ...
      'verbose',                 in.cvode_verbose);
    
  else
    [~, orbit_to_cycle_y] = feval(in.time_integration_method, ...
          @(t,y) dydt_ode(0, y, in.ode_parameters{:}), ...
          orbit_to_cycle_t, ...
          in.initial_point, ...
          in.time_integration_options);
  end
  
  point_on_cycle = orbit_to_cycle_y(end,:)';
  
  if in.show_plots
    my_figure = figure;
    plot_t = linspace(0, in.time_to_converge_to_cycle, ...
                         in.n_interpolated_points);
    transformed = in.plot_transformation(orbit_to_cycle_y);
    transformed(isinf(transformed)) = 0;
    transformed(isnan(transformed)) = 0;
    plot_y = interp1(orbit_to_cycle_t, transformed, plot_t, 'makima');
    plot(plot_t', plot_y);
    xlabel('t')
    ylabel(in.ylabel)
    disp('Now showing plot from  t = 0  to  t = time_to_converge_to_cycle')
    input('Press enter to continue or ctrl-c to abort', 's')
    if isvalid(my_figure)
      close(my_figure.Number)
    end
  end
end