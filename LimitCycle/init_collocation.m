
function initial_continuation_data = init_collocation(varargin)

  lc_input.point_on_limitcycle      = [];
  lc_input.odefile                  = [];
  lc_input.ode_parameters           = [];
  lc_input.active_parameter_index   = [];
  lc_input.time_integration_method  = @ode15s;
  lc_input.lower_bound_period       = [];
  lc_input.upper_bound_period       = [];
  lc_input.time_integration_options = odeset();
  lc_input.poincare_tolerance       = 5e-2;
  lc_input.show_plot                = false;
  lc_input.nCollocationPoints       = 4;
  lc_input.nMeshIntervals           = [];
  lc_input.collocation_tolerance    = 1e-2;
  
  i=1;
  while i <= nargin
    if ~ ischar(varargin{i})
      error('Please specify options as key-value pairs')
    end
    if ~ isfield(lc_input,varargin{i})
      error([varargin{i} ' is not a valid option.'])
    end
    lc_input.(varargin{i}) = varargin{i+1};
    i = i+2;
  end
  fields = fieldnames(lc_input);
  for i=1:length(fields)
    if isempty(lc_input.(fields{i}))
      error(['You must specifiy ' fields{i} '.'])
    end
  end

  [lc_input.t, lc_input.y] = compute_periodic_solution(lc_input);

  lc_input.y             = lc_input.y'; % transpose y
  lc_input.p             = cell2mat(lc_input.ode_parameters);
  lc_input.ap            = lc_input.active_parameter_index;
  lc_input.ntst          = lc_input.nMeshIntervals;
  lc_input.ncol          = lc_input.nCollocationPoints;
  lc_input.tolerance     = lc_input.collocation_tolerance;
  
  [initial_continuation_data, ~] = ...
    initOrbLC_L(lc_input.odefile, ...
                lc_input.t, ...
                lc_input.y, ...
                lc_input.p, ...      
                lc_input.ap, ...
                lc_input.ntst, ...
                lc_input.ncol, ...    
                lc_input.tolerance);
              
end

function [solution_t, solution_x] = compute_periodic_solution(in)
    
  global cds
  
  handles      = feval(in.odefile);
  dydt_ode     = handles{2};
  jacobian_ode = handles{3};
  cds.nphases  = length(in.point_on_limitcycle);

  tangent_to_limitcycle ...   
    = dydt_ode(0,in.point_on_limitcycle, in.ode_parameters{:});
  time_integration_options = odeset(in.time_integration_options, ...
    'Events',       @returnToPlane, ...
    'Jacobian',     @(t,y) jacobian_ode(0, y, in.ode_parameters{:}));
  
  solution = feval(in.time_integration_method, ...
    @(t,y) dydt_ode(0, y, in.ode_parameters{:}), ...
    [0,in.upper_bound_period], ...
    in.point_on_limitcycle, ...
    time_integration_options); 
  
  period                    = solution.x(end);
  
  if in.show_plot
    trail_solution_t = linspace(0, period, 500);
    trail_solution_x = deval(solution, trail_solution_t);
    plot(trail_solution_t, trail_solution_x-solution.y(:,1))
    xlabel('t')
    ylabel('deviation form initial value')
    pause
  end
  

  solution = odextend(solution,[],1.1*period);
  
  solution_t = linspace(0,1.1*period, ...
    3 * in.nMeshIntervals * in.nCollocationPoints);
  solution_x = deval(solution, solution_t);
  
  
  if in.show_plot
    plot(solution_t, solution_x - solution_x(:,1))
    xlabel('t')
    ylabel('deviation form initial value')
    pause
  end
  
  function [value, isterminal, direction] = returnToPlane(t, x)
    % x and should be a column vector
    value = tangent_to_limitcycle'*(x-in.point_on_limitcycle);
    isterminal = t > in.lower_bound_period ...
      && max(abs(x-in.point_on_limitcycle)) < in.poincare_tolerance;
    direction = 1;
  end
end
