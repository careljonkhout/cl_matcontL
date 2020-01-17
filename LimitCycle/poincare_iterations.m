function point = poincare_iterations(x)
  global cds contopts
  
  
  iteration                  = 1;
  lowest_curve_function_norm = Inf;
  done                       = false;
  while ~ done
    [point_on_limitcycle, period, ode_parameters] = ...
            NP_SS_extract_phases_period_and_parameters(x);
          
    tangent_to_limitcycle = ...   
            cds.dydt_ode(0, point_on_limitcycle, ode_parameters{:});
    
          
    abs_tol = contopts.integration_abs_tol;
    rel_tol = contopts.integration_rel_tol;
    
    end_time = contopts.PoincareMaxPeriodIncFactor * period;

    if cds.using_cvode
      [orbit_t, orbit_y] = cds.integrator( ...
              't_values',            [0, end_time], ...
              'initial_point',       point_on_limitcycle, ...
              'cycle_detection',     true, ...
              'lower_bound_period',  period / 10, ...
              'point_on_limitcycle', point_on_limitcycle, ...
              'ode_parameters',      cell2mat(ode_parameters), ...
              'abs_tol',             abs_tol, ...
              'rel_tol',             rel_tol);
    else
      time_integration_options = odeset( ...
              'AbsTol', abs_tol, ...
              'RelTol', rel_tol, ...
              'Events', @returnToPlane);
      if ~ isempty(cds.jacobian_ode)
        time_integration_options = odeset(in.time_integration_options, ...
                'Jacobian', @(t,y) cds.jacobian_ode(0, y, ode_parameters{:}));
      end
      [orbit_t, orbit_y] = cds.integrator( ...
        @(t,y) cds.dydt_ode(0, y, ode_parameters{:}), ...
        [0, end_time], ...
        point_on_limitcycle, ...
        time_integration_options);
    end
    
    cycle_detected = orbit_t(end) ~= end_time;
    
    if ~ cycle_detected
      print_diag(0, 'Poincare iterations failed. No cycle was detected.\n')
      break;
    end

    new_period              = orbit_t(end);
    new_point_on_limitcycle = orbit_y(end, :)';
    active_parameter_value  = x(end);
    
    new_x = [new_point_on_limitcycle; new_period; active_parameter_value];
    
    point_lc     = point_on_limitcycle;
    new_point_lc = new_point_on_limitcycle;
    
    curve_function_norm        = max(abs( new_point_lc - point_lc));
    lowest_curve_function_norm = min(lowest_curve_function_norm, ...
                                            curve_function_norm);
    
    

                                          
                                          
    rel_curve_function_norm = max(abs((new_point_lc - point_lc) ./ point_lc));
   
    print_diag(0, ['curve function norm: 10^%.2f ' ...
                   'relative curve function norm: 10^%.2f ', ...
                   'iteration: %d\n'], ...
                    log10(curve_function_norm), ...
                    log10(rel_curve_function_norm), ...
                    iteration);
                  
    if curve_function_norm > ...
               contopts.max_rel_funcnorm_increase * lowest_curve_function_norm
      print_diag(1,[ ...
        'Current curve function norm is now more than %d times ' ...
        'the lowest curve function norm that ' ...
        'was attained in this series of Poincare iterations. ' ...
        'Aborting Poincare Iterations.\n'], contopts.max_rel_funcnorm_increase);
      break;
    end
    done = curve_function_norm < contopts.FunTolerance || ...
           iteration           >= contopts.MaxCorrIters;
    iteration = iteration + 1;      
    x = new_x;
    

%     if in.show_plots
% 
%       my_figure = figure;
%       plot_t = linspace(0, period, in.n_interpolated_points);
%       transformed_orbit = in.plot_transformation(orbit_y) ...
%                         - in.plot_transformation(orbit_y(1,:));
%       plot_y = interp1(orbit_t, transformed_orbit, plot_t, in.interpolation);
%       plot(plot_t', plot_y);
%       xlabel('t')
%       ylabel('deviation form initial value')
%       disp(['Now showing plot from t=time_to_converge_to_cycle to ' ...
%                                          't=time_to_converge_to_cycle + period']);
%       my_pause();
%       if isvalid(my_figure)
%         close(my_figure.Number)
%       end
%     end

%     if ~ using_cvode
%       fprintf('max-norm of gap in cycle: %.4e\n', norm_of_gap)
%     end

  end
  
  if cycle_detected && curve_function_norm < contopts.FunTolerance
    point.R                            = curve_function_norm;
    point.relative_curve_function_norm = rel_curve_function_norm;
    point.x                            = x;
    point.iters                        = iteration;
  else
    point = [];
  end
 
  function [value, isterminal, direction] = returnToPlane(t, x)
    % x and should be a column vector
    value = tangent_to_limitcycle' * (x - point_on_limitcycle);
    isterminal = t > in.lower_bound_period;
    direction = 1;
  end

end