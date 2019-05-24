
function initial_continuation_data = init_collocation_find_stable_cycle(varargin)

  input.initial_point             = [];
  input.time_to_converge_to_cycle = [];
  input.odefile                   = [];
  input.ode_parameters            = [];
  input.active_parameter_index    = [];
  input.time_integration_method   = @ode15s;
  input.lower_bound_period        = [];
  input.upper_bound_period        = [];
  input.time_integration_options  = odeset();
  input.poincare_tolerance        = 5e-2;
  input.show_plot                 = false;
  input.nCollocationPoints        = 4;
  input.nMeshIntervals            = [];
  input.collocation_tolerance     = 1e-2;
  
  i=1;
  while i <= nargin
    if ~ ischar(varargin{i})
      error('Please specify options as name-value pairs')
    end
    if ~ isfield(input,varargin{i})
      error([varargin{i} ' is not a valid option.'])
    end
    input.(varargin{i}) = varargin{i+1};
    i = i+2;
  end
  fields = fieldnames(input);
  for i=1:length(fields)
    if isempty(input.(fields{i}))
      error(['You must specifiy ' fields{i} '.'])
    end
  end
  
  ncol = input.nCollocationPoints;
  
  if ~ isnumeric(ncol) || numel(ncol) ~= 1 || floor(ncol) ~= ncol || ...
       ncol < 2 || ncol > 7
    error('nCollocationPoints must be an integer between 1 and 8')
  end

  ntst = input.nMeshIntervals;
  
  if ~ isnumeric(ntst) || numel(ntst) ~= 1 || floor(ntst) ~= ntst || ntst < 2
    error('nMeshIntervals must be an integer greater than 1')
  end
  input.point_on_limitcycle = converge_to_cycle(input);
  %lc_input.point_on_limitcycle = lc_input.point_on_limitcycle';
  
  [input.t, input.y] = compute_periodic_solution(input);

  input.y             = input.y'; % transpose y
  input.p             = cell2mat(input.ode_parameters);
  input.ap            = input.active_parameter_index;
  input.ntst          = input.nMeshIntervals;
  input.ncol          = input.nCollocationPoints;
  input.tolerance     = input.collocation_tolerance;
  
  [initial_continuation_data, ~] = ...
    initOrbLC_L(input.odefile, ...
                input.t, ...
                input.y, ...
                input.p, ...      
                input.ap, ...
                input.ntst, ...
                input.ncol, ...    
                input.tolerance);
              
end

function point_on_cycle = converge_to_cycle(in)
    
  global cds
  
  handles      = feval(in.odefile);
  dydt_ode     = handles{2};
  jacobian_ode = handles{3};
  cds.nphases  = length(in.initial_point);

 
  if ~ isempty(jacobian_ode)
    in.time_integration_options = odeset(in.time_integration_options, ...
      'Jacobian',     @(t,y) jacobian_ode(0, y, in.ode_parameters{:}));
  end
  
  solution = feval(in.time_integration_method, ...
      @(t,y) dydt_ode(0, y, in.ode_parameters{:}), ...
      [0, in.time_to_converge_to_cycle], ...
      in.initial_point, ...
      in.time_integration_options); 
  
  if in.show_plot
    orbit_t = linspace(0, solution.x(end), 500);
    orbit_x = deval(solution, orbit_t);
    plot(orbit_t, orbit_x)
    plot(orbit_t, orbit_x-solution.y(:,1))
    xlabel('t')
    ylabel('phase variables')
    pause
  end
  
  point_on_cycle = deval(solution, in.time_to_converge_to_cycle);

end


function [solution_t, solution_x] = compute_periodic_solution(in)
    
  global cds
  
  handles      = feval(in.odefile);
  dydt_ode     = handles{2};
  jacobian_ode = handles{3};
  cds.nphases  = length(in.point_on_limitcycle);

  tangent_to_limitcycle ...   
    = dydt_ode(0,in.point_on_limitcycle, in.ode_parameters{:});
  in.time_integration_options = odeset(in.time_integration_options, ...
    'Events',       @returnToPlane);
  
  if ~ isempty(jacobian_ode)
    in.time_integration_options = odeset(in.time_integration_options, ...
      'Jacobian',     @(t,y) jacobian_ode(0, y, in.ode_parameters{:}));
  end
  
  solution = feval(in.time_integration_method, ...
    @(t,y) dydt_ode(0, y, in.ode_parameters{:}), ...
    [0,in.upper_bound_period], ...
    in.point_on_limitcycle, ...
    in.time_integration_options); 
  
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
