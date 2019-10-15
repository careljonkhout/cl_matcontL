function initial_continuation_data = init_EP_integration(varargin)

  contopts = contset();
  % required arguments
  input.initial_point                   = [];
  input.time_to_converge_to_equilibrium = [];
  input.odefile                         = [];
  input.ode_parameters                  = [];
  input.active_parameter_index          = [];
  
  % optional arguments, i.e. arguments with default values
  input.time_integration_method         = @ode15s;
  input.time_integration_options        = odeset( ...
    'AbsTol', contopts.integration_abs_tol, ...
    'RelTol', contopts.integration_rel_tol);
  input.show_plot                       = false;
  
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
  
  if ~ iscell(input.ode_parameters)
    input.ode_parameters = num2cell(input.ode_parameters);
  end
  
  equilibrium = converge_to_equilibrium(input);
  initial_continuation_data = init_EP_EP_L( ...
          input.odefile, ...
          equilibrium, ...
          input.ode_parameters, ...
          input.active_parameter_index); 
end
