function initial_continuation_data = init_multiple_shooting_from_pd(varargin)
                                                      
  input.odefile                   = [];
  input.period_doubling_point     = [];
  input.epsilon                   = [];
  input.time_integration_method   = @ode15s;
  input.basis_size                = [];
  input.show_plots                = false;
  
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
  nMeshIntervals        = length(pd.time_mesh) - 1;
  nphases               = (length(pd.x) - 2) / nMeshIntervals;
  point_on_limitcycle   = pd.x(1:nphases);
  ode_handles           = feval(input.odefile);
  dydt_ode              = ode_handles{2};
  using_cvode           = endsWith(func2str(time_integration_method),'cvode');
  basis_size            = input.basis_size;
  
  global cds;
  cds = [];
  cds.nphases                    = nphases;
  cds.probfile                   = odefile;
  cds.options.PartitionMonodromy = cds.nphases > 20;
  cds.nap                        = 1;
  cds.ndim                       = cds.nphases * nMeshIntervals + 2;
  cds.usernorm                   = [];
  cds.ncoo                       = cds.nphases * nMeshIntervals + 1;
  cds.ActiveParams               = pd.ap;
  cds.P0                         = cell2mat(parameters);
  cds.previous_phases            = point_on_limitcycle(:);
  cds.dydt_ode                   = dydt_ode;
  cds.jacobian_ode               = ode_handles{3};
  cds.jacobian_p_ode             = ode_handles{4};
  cds.integrator                 = time_integration_method;
  cds.preferred_basis_size       = basis_size;
  cds.p                          = basis_size;
  cds.mv_count                   = 0;
  cds.curve                      = @multiple_shooting;
  cds.nMeshIntervals             = nMeshIntervals;
  cds.mesh                       = pd.time_mesh;
  cds.using_cvode                = using_cvode;

  period    = pd.x(end-1);   
  pd_vector = NewtonPicard.MultipleShooting.find_pd_vector(pd.x);
  h         = input.epsilon;
  
  
  initial_continuation_data = init_multiple_shooting( ...
    'point_on_limitcycle',      pd.x(1:cds.nphases) + h * pd_vector, ...
    'odefile',                  input.odefile, ...
    'ode_parameters',           parameters, ...
    'active_parameter_index',   pd.ap, ...
    'lower_bound_period',       1.9 * period, ...
    'upper_bound_period',       2.1 * period, ...
    'nMeshIntervals',           pd.nMeshIntervals, ...
    'subspace_size',            input.basis_size, ...
    'time_integration_method',  input.time_integration_method, ...
    'show_plots',               input.show_plots);
  
end