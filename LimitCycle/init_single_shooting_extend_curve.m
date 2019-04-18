
function init_single_shooting_extend_curve(varargin)

  lc_input.initial_continuation_state          = [];
  lc_input.initial_continuation_tangent_vector = [];
  lc_input.odefile                             = [];
  lc_input.ode_parameters                      = [];
  lc_input.active_parameter_index              = [];
  lc_input.time_integration_method             = @ode15s;
  lc_input.time_integration_options            = odeset();
  lc_input.subspace_size                       = [];
  lc_input.nDiscretizationPoints               = 100;
  
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

  
  do_init_single_shooting_extend_curve(lc_input);
              
end

function do_init_single_shooting_extend_curve(in)
    
  global cds  
  handles                = feval(in.odefile);
  dydt_ode               = handles{2};
  jacobian_ode           = handles{3};
  cds.nphases            = length(in.initial_continuation_state) - 2;
  
  point_on_limitcycle    = in.initial_continuation_state(1:end-2);
  tangent_to_limitcycle  = dydt_ode(0,point_on_limitcycle,in.ode_parameters{:});

  cds.probfile        = in.odefile;
  cds.options.PartitionMonodromy = cds.nphases > 20;
  cds.nap             = 1;
  cds.ndim            = cds.nphases + 2;
  cds.usernorm        = [];
  cds.ncoo            = cds.nphases + 1;
  cds.ActiveParams    = in.active_parameter_index;
  cds.P0              = cell2mat(in.ode_parameters);
  cds.previous_phases = point_on_limitcycle;
  cds.previous_dydt_0 = tangent_to_limitcycle;
  cds.dydt_ode        = dydt_ode;
  cds.jacobian_ode    = jacobian_ode;
  cds.integrator      = in.time_integration_method;
  cds.preferred_basis_size  = in.subspace_size;
  cds.nDiscretizationPoints = in.nDiscretizationPoints;
  cds.p               = in.subspace_size;
  cds.mv_count        = 0;
 
end
