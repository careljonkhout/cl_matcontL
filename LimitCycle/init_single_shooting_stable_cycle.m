
function initial_continuation_data = init_single_shooting_stable_cycle(varargin)

  lc_input.initial_point             = [];
  lc_input.time_to_converge_to_cycle = [];
  lc_input.odefile                   = [];
  lc_input.ode_parameters            = [];
  lc_input.active_parameter_index    = [];
  lc_input.time_integration_method   = @ode15s;
  lc_input.lower_bound_period        = [];
  lc_input.upper_bound_period        = [];
  lc_input.time_integration_options  = odeset();
  lc_input.poincare_tolerance        = 5e-2;
  lc_input.show_plot                 = false;
  lc_input.subspace_size             = [];
  lc_input.nDiscretizationPoints     = 100;
  
  i=1;
  while i <= nargin
    if ~ ischar(varargin{i})
      error('Please specify options as name-value pairs')
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

  lc_input.point_on_limitcycle = converge_to_cycle(lc_input);
  %lc_input.point_on_limitcycle = lc_input.point_on_limitcycle';
  
  initial_continuation_data = do_init_single_shooting(lc_input);
              
end

function point_on_cycle = converge_to_cycle(in)
    
  global cds
  
  handles      = feval(in.odefile);
  dydt_ode     = handles{2};
  jacobian_ode = handles{3};
  cds.nphases  = length(in.initial_point);

 
  time_integration_options = odeset(in.time_integration_options, ...
    'Jacobian',     @(t,y) jacobian_ode(0, y, in.ode_parameters{:}));
  
  solution = feval(in.time_integration_method, ...
      @(t,y) dydt_ode(0, y, in.ode_parameters{:}), ...
      [0, in.time_to_converge_to_cycle], ...
      in.initial_point, ...
      time_integration_options); 
  
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

