function point_on_cycle = converge_to_cycle(in)
    
  global cds;
  
  handles      = feval(in.odefile);
  dydt_ode     = handles{2};
  jacobian_ode = handles{3};
  cds.nphases  = length(in.initial_point);
  using_cvode  = endsWith(func2str(in.time_integration_method), 'cvode');
  
  if ~ isempty(jacobian_ode) && ~ using_cvode
    in.time_integration_options = odeset(in.time_integration_options, ...
                   'Jacobian', @(t,y) jacobian_ode(0, y, in.ode_parameters{:}));
  end
  
  orbit_to_cycle_t = linspace(0, in.time_to_converge_to_cycle, 5000);
  
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
    plot(orbit_to_cycle_t, orbit_to_cycle_y);
    xlabel('t')
    ylabel('phase variables')
    disp('Now showing plot from t = 0 to t = time_to_converge_to_cycle')
    disp('Press a key to continue')
    pause
    if isvalid(my_figure)
      close(my_figure.Number)
    end
  end
end