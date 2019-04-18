
function initial_continuation_data = ...
    init_single_shooting_find_stable_cycle(varargin)
  
  contopts = contset();

  input.initial_point             = [];
  input.time_to_converge_to_cycle = [];
  input.odefile                   = [];
  input.ode_parameters            = [];
  input.active_parameter_index    = [];
  input.time_integration_method   = @ode15s;
  input.lower_bound_period        = [];
  input.upper_bound_period        = [];
  input.time_integration_options = odeset( ...
    'AbsTol', contopts.integration_abs_tol, ...
    'RelTol', contopts.integration_rel_tol);
  input.poincare_tolerance        = 5e-2;
  input.show_plot                 = false;
  input.subspace_size             = [];
  input.nDiscretizationPoints     = 100;
  
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

  input.point_on_limitcycle = converge_to_cycle(input);
  
  initial_continuation_data = init_single_shooting_internal(input);    
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
    orbit_to_cycle_t = linspace(0, in.time_to_converge_to_cycle, 500);
    orbit_to_cycle_x = deval(solution, orbit_to_cycle_t);
    plot(orbit_to_cycle_t, orbit_to_cycle_x-solution.y(:,1))
    xlabel('t')
    ylabel('phase variables')
    pause
  end
  
  point_on_cycle = deval(solution, in.time_to_converge_to_cycle);

end

