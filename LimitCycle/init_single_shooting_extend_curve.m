function init_single_shooting_extend_curve(varargin)

  input.initial_continuation_state          = [];
  input.odefile                             = [];
  input.ode_parameters                      = [];
  input.active_parameter_index              = [];
  input.time_integration_method             = @ode15s;
  input.subspace_size                       = [];
  
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

  
  do_init_single_shooting_extend_curve(input);
              
end

function do_init_single_shooting_extend_curve(in)
  
  clear global
  global cds  
  handles                = feval(in.odefile);
  dydt_ode               = handles{2};
  jacobian_ode           = handles{3};
  cds.n_phases           = length(in.initial_continuation_state) - 2;
  ode_parameters         = convert_to_cell_if_needed(in.ode_parameters);
  
  point_on_limitcycle    = in.initial_continuation_state(1:end-2);
  tangent_to_limitcycle  = dydt_ode(0,point_on_limitcycle,ode_parameters{:});
  
  
  cds.using_cvode     = false; % todo: add cvode support
  cds.probfile        = in.odefile;
  cds.options.PartitionMonodromy = cds.n_phases > 20;
  cds.nap             = 1;
  cds.ndim            = cds.n_phases + 2;
  cds.usernorm        = [];
  cds.ncoo            = cds.n_phases + 1;
  cds.ActiveParams    = in.active_parameter_index;
  cds.P0              = cell2mat(ode_parameters);
  cds.previous_phases = point_on_limitcycle;
  cds.previous_dydt_0 = tangent_to_limitcycle;
  cds.dydt_ode        = dydt_ode;
  cds.jacobian_ode    = jacobian_ode;
  cds.jacobian_p_ode  = handles{4};
  cds.integrator      = in.time_integration_method;
  cds.preferred_basis_size  = in.subspace_size;
  cds.p               = in.subspace_size;
  cds.mv_count        = 0;
  cds.curve           = @single_shooting;
 
end
