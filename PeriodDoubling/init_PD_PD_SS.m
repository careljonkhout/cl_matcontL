function x = init_PD_PD_SS(varargin)
                                                      
  input.odefile                   = [];
  input.period_doubling_point     = [];
  input.epsilon                   = [];
  input.time_integration_method   = @ode15s;
  input.basis_size                = [];
  input.active_parameter_indices  = [];
  
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
  
  % todo: if ~ isfield(input.period_doubling_point, .....  raise error
  
  pd                    = input.period_doubling_point;
  odefile               = input.odefile;
  time_integration_method = input.time_integration_method;
  parameters            = pd.P0;
  parameters(pd.ap)     = pd.x(end);
  parameters            = num2cell(parameters);
  point_on_limitcycle   = pd.x(1:end-2);
  ode_handles           = feval(input.odefile);
  dydt_ode              = ode_handles{2};
  tangent_to_limitcycle = feval(dydt_ode,0, point_on_limitcycle, parameters{:});
  using_cvode           = endsWith(func2str(time_integration_method), 'cvode');
  basis_size            = input.basis_size;
  
  global cds;
  cds = [];
  cds.n_phases                    = length(pd.x) - 2;
  cds.probfile                   = odefile;
  cds.options.PartitionMonodromy = cds.n_phases > 20;
  cds.nap                        = 2;
  cds.ndim                       = 2 * cds.n_phases + 3;
  cds.usernorm                   = [];
  cds.ncoo                       = 2 * cds.n_phases + 1;
  cds.P0                         = cell2mat(parameters);
  cds.ActiveParams               = pd.ap;
  cds.previous_phases            = point_on_limitcycle(:);
  cds.previous_dydt_0            = tangent_to_limitcycle(:);
  cds.dydt_ode                   = dydt_ode;
  cds.jacobian_ode               = ode_handles{3};
  cds.jacobian_p_ode             = ode_handles{4};
  cds.integrator                 = time_integration_method;
  cds.preferred_basis_size       = basis_size;
  cds.p                          = basis_size;
  cds.mv_count                   = 0;
  cds.curve                      = @perioddoubling_ss;
  cds.using_cvode                = using_cvode;
  period                         = pd.x(end-1);
  
  
  pd_vector = NP_SS_find_pd_vector(pd.x);
  
  cds.l = pd_vector;
  
  cds.ActiveParams               = input.active_parameter_indices;
  active_parameter_values = ...
                           cell2mat(parameters(input.active_parameter_indices));
  
  x = [pd.x(1:end-2); pd_vector; period; active_parameter_values(:)];
end