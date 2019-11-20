%% initialize a cycle continuation with collocation from a stable cycle.
% This function will handle the integration towards the stable cycle. The
% arguments to this function must be specified as name value pairs, for
% instance:
%
% initial_continuation_data = init_collocation_find_stable_cycle( ...
%   'initial_point',             ones(100,1), ...
%   'time_to_converge_to_cycle', 150, ... 
%   'odefile',                   @my_odefile,  ...
%   'ode_parameters',            {0.1, 0.2, 0.3}, ...
%   'active_parameter_index',    3, ...
%   'lower_bound_period',        1, ...
%   'upper_bound_period',        20, ...
%   'n_mesh_intervals',            20, ...
%   'show_plots',                false ...
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
%% subspace_size
% the preferred size of the initial subspace used in the Newton-Picard method.
% Note that the subspace size in continuously adjusted. See contset.
%
%% +++++ Optional arguments +++++
%
%% time_integration_method
% The matlab odesuite time integration function, or a compatible alternative
% that is to be used in the finding of the limit cycle initialization.
% Collocation does not use integration, so this setting does not influence the
% continuantion, only the initialization. Default value: @ode15s
%
%% time_integration_options
% the time integration options that are to be used during the finding of the
% limit cycle. Collocation does not use integration, so this settings do not
% influence the continuantion, only the initialization. Default value:
% odeset('AbsTol', contopts.integration_abs_tol, ...
%        'RelTol', contopts.integration_rel_tol);
% where contopts contain the default options from contset.
%
%% poincare_tolerance
% The maximum distance between to consecutive iterates of the Poincare map that
% is still considered to be a cycle. See the code in
% init_multiple_shooting_internal and "how a cycle is detected" above to see how
% poincare_tolerance is used internally. Default value: 2e-3
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
%
%% collocation_tolerance
% This function calls the function initOrbLC_L, which also has a tolerance
% parameter "h" related to the gap in the almost converged cycle, similar to
% poincare_tolerance. The value of collocation_tolerance is passed on to
% initOrbLC as "h".Default value: 1e-2
%
%% nCollocationPoints
% the degree of the polynomials that are used to represent the cycle. This
% parameter must be at least 2 and at most 7. The default value of 4 seems to
% always work well. The number of nonzero's in the Jacobian matrix of the
% continuation increases as the square of nCollocationPoints, so higher values
% of nCollocationPoints will likely increase computation times significantly.
% Therefore, if one wants to increase the accuracy of the cycles in the
% continuation, it is recommended to increase n_mesh_intervals to increase the
% accuracy of the continuation, and to keep nCollocationPoints at the default
% value (, although it has not been ruled out entirely that higher values of
% nCollocationPoints might be advantageous in some special cases.) 
% Default value: 4

function initial_continuation_data = init_collocation_find_stable_cycle(varargin)

  contopts = contset();
  % required arguments
  input.initial_point             = [];
  input.time_to_converge_to_cycle = [];
  input.odefile                   = [];
  input.ode_parameters            = [];
  input.active_parameter_index    = [];
  input.lower_bound_period        = [];
  input.upper_bound_period        = [];
  input.n_mesh_intervals          = [];
  input.plot_coordinates          = 'all';
  input.cvode_verbose             = false;
  
  % optional arguments, i.e. arguments with default values
  input.time_integration_method   = @ode15s;
  input.time_integration_options = odeset( ...
    'AbsTol', contopts.integration_abs_tol, ...
    'RelTol', contopts.integration_rel_tol);
  input.poincare_tolerance        = 2e-3;
  input.show_plots                = false;
  input.collocation_tolerance     = 1e-2;
  input.nCollocationPoints        = 4;
  input.plot_transformation       = @(x) x;
  input.ylabel                    = 'phase variables';
  input.n_computed_points         = 100;
  input.n_interpolated_points     = 10000;
  input.interpolation             = 'makima';
  
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
  
  ncol = input.nCollocationPoints;
  
  if ~ isnumeric(ncol) || numel(ncol) ~= 1 || floor(ncol) ~= ncol || ...
       ncol < 2 || ncol > 7
    error('nCollocationPoints must be an integer, at least 2, and at most 7')
  end

  ntst = input.n_mesh_intervals;
  
  if ~ isnumeric(ntst) || numel(ntst) ~= 1 || floor(ntst) ~= ntst || ntst < 2
    error('n_mesh_intervals must be an integer greater than 1')
  end
  
  if ~ iscell(input.ode_parameters)
    input.ode_parameters = num2cell(input.ode_parameters);
  end
  
  
  input.point_on_limitcycle = converge_to_cycle(input);
 
 
  
  [input.t, input.y] = compute_periodic_solution(input);

  input.y             = input.y'; % transpose y
  input.p             = cell2mat(input.ode_parameters);
  input.ap            = input.active_parameter_index;
  input.ntst          = input.n_mesh_intervals;
  input.ncol          = input.nCollocationPoints;
  input.tolerance     = input.collocation_tolerance;
  
  [initial_continuation_data, ~] = ...
    initOrbLC_L(input.odefile, ...
                input.t, ...
                input.y, ...
                input.p, ...      
                input.ap, ...
                input.ntst, ...
                input.ncol, ...    
                input.tolerance);
              
end

function [solution_t, solution_x] = compute_periodic_solution(in)
    
  global cds
  
  handles      = feval(in.odefile);
  dydt_ode     = handles{2};
  jacobian_ode = handles{3};
  cds.n_phases  = length(in.point_on_limitcycle);

  tangent_to_limitcycle = dydt_ode( ...
                               0, in.point_on_limitcycle, in.ode_parameters{:});
                     
  in.time_integration_options = odeset(in.time_integration_options, ...
                                        'Events',       @returnToPlane);
  
  if ~ isempty(jacobian_ode)
    in.time_integration_options = odeset(in.time_integration_options, ...
      'Jacobian',     @(t,y) jacobian_ode(0, y, in.ode_parameters{:}));
  end
  
  solution = feval(in.time_integration_method, ...
    @(t,y) dydt_ode(0, y, in.ode_parameters{:}), ...
    [0,in.upper_bound_period], ...
    in.point_on_limitcycle, ...
    in.time_integration_options); 
  
  period   = solution.x(end);
  solution = odextend(solution, [], 1.1 * period);
  
  solution_t = linspace(0, 1.1 * period, ...
          3 * in.n_mesh_intervals * in.nCollocationPoints);
  solution_x = deval(solution, solution_t);
  
  
  if in.show_plots
    my_figure = figure;
    if strcmp(in.plot_coordinates,'all')
      plot(solution_t, solution_x - solution_x(:,1))
    else
      plot(solution_t, solution_x(in.plot_coordinates, :) ...
                     - solution_x(in.plot_coordinates, 1))
    end
    xlabel('t')
    ylabel('deviation form initial value')
    disp(['Now showing plot from t = time_to_converge_to_cycle to ' ...
                               't = time_to_converge_to_cycle + 1.1 * period']);
    my_pause();
    if isvalid(my_figure)
      close(my_figure.Number)
    end
  end
  
  function [value, isterminal, direction] = returnToPlane(t, x)
    % x and should be a column vector
    value = tangent_to_limitcycle'*(x-in.point_on_limitcycle);
    isterminal = t > in.lower_bound_period ...
      && max(abs(x-in.point_on_limitcycle)) < in.poincare_tolerance;
    direction = 1;
  end
end
