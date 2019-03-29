%
% Initializes the global struct variable cds, and computes
% initial_continuation_data from an orbit that has converged to a limit cycle
% supplied by the user.
% 
% values must be specified as key, value pairs, for instance:
%
% initial_continuation_data = init_single_shooting( ...
%  'point_on_limitcycle',     < point on limit cycle >, ...
%  'odefile',                 < odefile >,  ...
%  'ode_parameters',          < cell array of ode parameters >, ...
%  'active_parameter_index',  < active parameter index >, ...
%  'lower_bound_period',      < lower bound period >, ...
%  'upper_bound_period',      < upper bound period >, ...
%  'subspace_size',           6, ...
%  'nDiscretizationPoints',   1000, ...
%  'show_plot'            ,   false);
% 
% where you should replace the things between < and > be the appropriate data
% from your problem
%
% inputs:
% 
% point_on_limitcycle -- a point on the limitcycle
%
% odefile -- The odefile with the definition of the system of ODE's in standard
% matcontL format.
% 
% ode_parameters -- a cell array of the values of the parameters of the system
% of ODE's for which the orbit that converges to the limitcycle was found
%
% active_parameter_index -- the index ( in ode_parameters ) of the parameter
% that is varied in the continuation.
%
% time_integration_method -- one of the methods of the matlab ode suite, or a
% fully compatible alternative. (ode15s seems to works well for spatially
% discretized pde's)
%
% lower_bound_period -- a lower bound for the period of the cycle.
%
% upper_bound_period -- an upperbound for the period of the cycle.
%
% time_integration_options -- the options structure that are supplied to the
% time_integation_method. Use the matlab command odeset to generate this
% structure. See the documentation of odeset.
%
% poincare_tolerance -- The cycle is determined to be "closed", if the
% max-norm-distance from the starting point is less than poincare_tolerance (the
% end point of orbit). If one is trying to start a cycle that loops around twice
% before closing, one might want to set this to a small value, and supply very
% small tolerances to time_integration_options. Otherwise, don't set it too
% small (e.g. not smaller than 1e-2 for cycle with amplitudes between 1 and 10),
% or the closing of the cycle might not be detected.



% allows for default values for init_single_shooting to be set in a nice way
function initial_continuation_data = init_single_shooting(varargin)

  ss_input.point_on_limitcycle      = [];
  ss_input.odefile                  = [];
  ss_input.ode_parameters           = [];
  ss_input.active_parameter_index   = [];
  ss_input.time_integration_method  = @ode15s;
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
    ss_input.(varargin{i}) = varargin{i+1};
    i = i+2;
  end
  fields = fieldnames(ss_input);
  for i=1:length(fields)
    if isempty(ss_input.(fields{i}))
      error(['You must specifiy ' fields{i} '.'])
    end
  end
  initial_continuation_data = do_init_single_shooting(ss_input);
end


function initial_continuation_data = do_init_single_shooting(in)
    
  global cds
  
  handles      = feval(in.odefile);
  dydt_ode     = handles{2};
  jacobian_ode = handles{3};
  cds.nphases  = length(in.point_on_limitcycle);

  tangent_to_limitcycle ...   
    = dydt_ode(0,in.point_on_limitcycle, in.ode_parameters{:});
  time_integration_options = odeset(in.time_integration_options, ...
    'Events',       @returnToPlane, ...
    'Jacobian',     @(t,y) jacobian_ode(0, y, in.ode_parameters{:}));
  
  [orbit_t, orbit_x] = feval(in.time_integration_method, ...
    @(t,y) dydt_ode(0, y, in.ode_parameters{:}), ...
    linspace(0,in.upper_bound_period), ...
    in.point_on_limitcycle, ...
    time_integration_options); 
  
  if in.show_plot
    plot(orbit_t, orbit_x-orbit_x(1,:))
    xlabel('t')
    ylabel('deviation form initial value')
    pause
  end
  
  period                    = orbit_t(end);
  initial_continuation_data = zeros(cds.nphases + 2, 1);
  initial_continuation_data(1:cds.nphases) = in.point_on_limitcycle;
  
  
  initial_continuation_data(end-1) = period;
  initial_continuation_data(end) = in.ode_parameters{in.active_parameter_index};
  
  cds.probfile        = in.odefile;
  cds.options.PartitionMonodromy = cds.nphases > 20;
  cds.nap             = 1;
  cds.ndim            = cds.nphases + 2;
  cds.usernorm        = [];
  cds.ncoo            = cds.nphases + 1;
  cds.ActiveParams    = in.active_parameter_index;
  cds.P0              = cell2mat(in.ode_parameters);
  cds.previous_phases = in.point_on_limitcycle;
  cds.previous_dydt_0 = tangent_to_limitcycle;
  cds.dydt_ode        = dydt_ode;
  cds.jacobian_ode    = jacobian_ode;
  cds.integrator      = in.time_integration_method;
  cds.preferred_basis_size  = in.subspace_size;
  cds.nDiscretizationPoints = in.nDiscretizationPoints;
  cds.p               = in.subspace_size;
  cds.mv_count        = 0;
 
  function [value, isterminal, direction] = returnToPlane(t, x)
    % x and should be a column vector
    value = tangent_to_limitcycle'*(x-in.point_on_limitcycle);
    isterminal = t > in.lower_bound_period ...
      && max(abs(x-in.point_on_limitcycle)) < in.poincare_tolerance;
    direction = 1;
  end
end