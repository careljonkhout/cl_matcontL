% allows for default values for init_single_shooting to be set in a nice way
function ss_input = sing_shoot_init_set(varargin)

  ss_input.point_on_limitcycle      = [];
  ss_input.odefile                  = [];
  ss_input.ode_parameters           = [];
  ss_input.active_parameter_index   = []; ...
  ss_input.time_integration_method  = @ode15s; ...
  ss_input.lower_bound_period       = [];
  ss_input.upper_bound_period       = [];
  ss_input.time_integration_options = odeset();
  ss_input.poincare_tolerance       = 1e-2;
  ss_input.nDiscretizationPoints    = 100;
  ss_input.subspace_size            = [];
  ss_input.show_plot                = false;
  
  i=1;
  while i <= nargin
    if ~ ischar(varargin{i})
      error('Please specify options as key-value pairs')
    end
    if ~ isfield(ss_input,varargin{i})
      error([varargin{i} ' is not a valid option.'])
    end
    eval(['ss_input.' varargin{i} '=varargin{i+1}'])
    i = i+2;
  end
  fields = fieldnames(ss_input);
  for i=1:length(fields)
    if isempty(getfield(ss_input,fields{i})) %#ok<GFLD>
      error(['You must specifiy ' fields{i} '.'])
    end
  end
  init
end