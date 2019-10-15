function x = init_single_shooting_from_pd(varargin)
                                                      
  input.odefile                   = [];
  input.period_doubling_point     = [];
  input.epsilon                   = [];
  input.time_integration_method   = @ode15s;
  input.basis_size                = [];
  
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
  h                     = input.epsilon;
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
  cds.nphases                    = length(pd.x) - 2;
  cds.probfile                   = odefile;
  cds.options.PartitionMonodromy = cds.nphases > 20;
  cds.nap                        = 1;
  cds.ndim                       = cds.nphases + 2;
  cds.usernorm                   = [];
  cds.ncoo                       = cds.nphases + 1;
  cds.ActiveParams               = pd.ap;
  cds.P0                         = cell2mat(parameters);
  cds.previous_phases            = point_on_limitcycle(:);
  cds.previous_dydt_0            = tangent_to_limitcycle(:);
  cds.dydt_ode                   = dydt_ode;
  cds.jacobian_ode               = ode_handles{3};
  cds.jacobian_p_ode             = ode_handles{4};
  cds.integrator                 = time_integration_method;
  cds.preferred_basis_size       = basis_size;
  cds.p                          = basis_size;
  cds.mv_count                   = 0;
  cds.curve                      = @single_shooting;
  cds.using_cvode                = using_cvode;
  
  error('this function is not yet completely implemented');
  
  period = pd.x(end-1);
  
  init_multiple_shooting( ...
    point_on_limitcycle          = pd.x(1:cds.nphases);
    input.odefile                = in.odefile;
    input.ode_parameters         = in.ode_parameters;
    input.active_parameter_index = pd.ap;
    input.lower_bound_period     = 1.9 * period;
    input.upper_bound_period     = 2.1 * period;
    input.nMeshIntervals         = pd.nMeshIntervals;
    input.subspace_size          = input.basis_size; 
 
    input.time_integration_method  = @ode15s;
    input.time_integration_options = odeset( ...
      'AbsTol', contopts.integration_abs_tol, ...
      'RelTol', contopts.integration_rel_tol);
    input.poincare_tolerance       = 1e-2;
    input.show_plots               = false;
    
    x_doubled = [x(1:end-2) + h * pd_vector; 2 * period; x(end)];
   x = NewtonPicard.SingleShooting.find_doubled_cycle(pd.x, h);
  
  
end