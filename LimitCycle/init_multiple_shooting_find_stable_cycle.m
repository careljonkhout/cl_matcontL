%% initialize a cycle continuation with multiple shooting from a stable cycle.
% This function will handle the integration towards the stable cycle. The
% arguments to this function must be specified as name value pairs, for
% instance:
%
% initial_continuation_data = init_multiple_shooting_find_stable_cycle( ...
%   'initial_point',             ones(100,1), ...
%   'time_to_converge_to_cycle', 150, ... 
%   'odefile',                   @my_odefile,  ...
%   'ode_parameters',            {0.1, 0.2, 0.3}, ...
%   'active_parameter_index',    3, ...
%   'n_mesh_intervals',            100, ...
%   'lower_bound_period',        1, ...
%   'upper_bound_period',        20, ...
%   'subspace_size',             10, ...
%   'show_plots',                 false ...
% );

%
%% +++++ how a cycle is detected +++++
% The integrator first integrates for a time of time_to_converge_to_cycle to the
% point p0. The Poincare plane is then defined as the plane perpendicular p0.
% Then the integration is continued until a return ( in the same direction ) to
% this Poincare plane is detected (for a maximum time interval of
% upper_bound_period ). Denote this second point of intersection of the
% trajectory with the Poincare plane by p1. If the max_norm of the difference
% between p0 and p1 is less than poincare_tolerance and the time from p0 to p1
% is larger than lower_bound_period, then the integration is stopped, and the
% period of the cycle is appoximated as the time it takes to go from p0 to p1.
% Default value: 1e-2
%
%% +++++ required arguments +++++
%
%% initial_point
% The initial point from which to start the time integration towards the stable
% cycle. 
%
%% time_to_converge_to_cycle
% The time is takes for the time integration to
% converge to the stable cycle enough, such that the difference between to
% subsequent iterates of the Poincare map is less than poincare_tolerance. The
% default value of poincare_tolerance is 1e-2.
%
%% odefile
% a function handle of the odefile that specifies the system of ODEs. For
% multiple shooting, it is currently required that the odefile must specify the
% Jacobian of the system of ODEs. If the Jacobian of the system of ODEs is not
% available, one could try to use collocation instead.
%
%% ode_parameters
% the parameter values of the system of ODEs, at which a stable cycle in present
%
%% active_parameter_index
% the 1-based index of the parameter that is to be varied during the
% continuation
%
%% lower_bound_period
% A lower bound on the period of the cycle. Setting this to high will result in
% two periods of the cycle being selected for the cycle. This is not desirable.
% See the code in init_multiple_shooting_internal and "how a cycle is detected"
% above to see how lower_bound_period is used internally.
%
%% upper_bound_period
% An upper bound on the period of the cycle. This is not desirable. See the code
% in init_multiple_shooting_internal and "how a cycle is detected" above to see
% how lower_bound_period is used internally.
%
%% n_mesh_intervals
% the number of mesh intervals used in the multiple shooting method.
%
%% subspace_size
% the preferred size of the initial subspace used in the Newton-Picard method.
% Note that the subspace size in continuously adjusted. See contset.
%
%% +++++ Optional arguments +++++
%
%% time_integration_method
% The matlab odesuite time integration function, or a compatible alternative
% that is to be used in the finding of the limit cycle initialization and the
% continuation. Default value: @ode15s
%
%% time_integration_options
% the time integration options that are to be used during the finding of the
% limit cycle. ( These options are NOT used during continuation. To set
% integration options for continuation, see contset). Default value:  
% odeset('AbsTol', contopts.integration_abs_tol, ...
%        'RelTol', contopts.integration_rel_tol);
% contopts contain the default options from contset.
%
%% poincare_tolerance
% The maximum distance between to consecutive iterates of the Poincare map that
% is still considered to be a cycle. See the code in
% init_multiple_shooting_internal and "how a cycle is detected" above to see how
% poincare_tolerance is used internally. Default value: 1e-2
%
%% show_plots
% set to true to see a plots of the ODE coordinates versus time. Two plots will
% be shown. The first plot starts form initial_points and has a time span of
% time_to_converge_to_cycle. The second plot has a time span between
% lower_bound_period and upper_bound_period, and start form the end point of the
% first plot. The second plot show the deviation from each variable from the
% first point of this plot, so that one can see the time point where the cycle
% return to its starting point as the point where all graph lines cross zero.
% After each plot is drawn the program pauses until a key is pressed. If
% show_plots is false, not plots are shown, and the program will not pause.
% Default value: false.

function initial_continuation_data = ...
          init_multiple_shooting_find_stable_cycle(varargin)
        
	contopts = contset();

  input.initial_point             = [];
  input.time_to_converge_to_cycle = [];
  input.odefile                   = [];
  input.ode_parameters            = [];
  input.active_parameter_index    = [];
  input.lower_bound_period        = [];
  input.upper_bound_period        = [];
  input.n_mesh_intervals            = [];
  input.subspace_size             = []; 
  input.time_integration_method   = @ode15s;
  input.time_integration_options  = odeset( ...
    'AbsTol', contopts.integration_abs_tol, ...
    'RelTol', contopts.integration_rel_tol);
  input.poincare_tolerance        = 1e-2;
  input.show_plots                = false;
  input.cvode_verbose             = false;
  input.plot_transformation       = @(x) x;
  input.ylabel                    = 'phase variables';
  input.n_computed_points         = 100;
  input.n_interpolated_points     = 10000;
  
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
  
  try
  Assert.scalar  ('time_to_converge_to_cycle', input.time_to_converge_to_cycle);
  Assert.positive('time_to_converge_to_cycle', input.time_to_converge_to_cycle);
  
  Assert.function_handle('odefile', input.odefile);
  
  Assert.scalar  ('active_parameter_index', input.active_parameter_index);
  Assert.positive('active_parameter_index', input.active_parameter_index);
  Assert.integer ('active_parameter_index', input.active_parameter_index);
  
  Assert.scalar  ('lower_bound_period', input.time_to_converge_to_cycle);
  Assert.positive('lower_bound_period', input.time_to_converge_to_cycle);
  
  Assert.scalar  ('upper_bound_period', input.time_to_converge_to_cycle);
  Assert.positive('upper_bound_period', input.time_to_converge_to_cycle);
  
  Assert.scalar      (   'n_mesh_intervals', input.n_mesh_intervals); 
  Assert.integer     (   'n_mesh_intervals', input.n_mesh_intervals);
  Assert.greater_than(1, 'n_mesh_intervals', input.n_mesh_intervals);

  Assert.scalar   ('subspace_size', input.subspace_size); 
  Assert.integer  ('subspace_size', input.subspace_size);
  Assert.positive ('subspace_size', input.subspace_size);
  assert(input.subspace_size <= length(input.initial_point), ...
      ['the subspace size must be less than or equal to ' ...
       'the dimension of the initial point']);
  catch failed_assertion
    error(failed_assertion.message);
  end
  
  if ~ iscell(input.ode_parameters)
    input.ode_parameters = num2cell(input.ode_parameters);
  end
    
  input.point_on_limitcycle = converge_to_cycle(input);
 
  initial_continuation_data = init_multiple_shooting_internal(input);             
end