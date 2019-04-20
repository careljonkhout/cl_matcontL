%
% Initializes the global struct variable cds, and computes
% initial_continuation_data from an orbit that has converged to a limit cycle
% supplied by the user.
% 
% inputs:
% 
% orbit -- An orbit suplied by the user that has converged to a limit cycle.
% Only the last point of the orbit is used. The second dimension of the orbit
% array should correspond the dimensions of the system of ode's.
%
% in.nMeshPoints -- The number of points used in multiple shooting 
%
% in.odefile -- The in.odefile with the definition of the system of ODE's in standard
% matcontL format.
% 
% in.ode_parameters -- a cell array of the values of the parameters of the system
% of ODE's for which the orbit that converges to the limitcycle was found
%
% in.active_parameter_index -- the index ( in in.ode_parameters ) of the parameter
% that is varied in the continuation.
%
% in.time_integration_method -- one of the methods of the matlab ode suite, or a
% fully compatible alternative. (ode15s seems to works well for spatially
% discretized pde's)
%
% in.lower_bound_period -- a lower bound for the period of the cycle.
%
% in.upper_bound_period -- an upperbound for the period of the cycle.
%
% in.time_integration_options -- the options structure that are supplied to the
% time_integation_method. Use the matlab command odeset to generate this
% structure. See the documentation of odeset.
%
% in.poincare_tolerance -- The cycle is determined to be "closed", if the
% max-norm-distance from the starting point is less than poincare_tolerance (the
% end point of orbit). If one is trying to start a cycle that loops around twice
% before closing, one might want to set this to a small value, and supply very
% small tolerances to in.time_integration_options. Otherwise, don't set it too
% small (e.g. not smaller than 1e-2 for cycle with amplitudes between 1 and 10),
% or the closing of the cycle might not be detected.
function initial_continuation_data = ...
        init_multiple_shooting(varargin)
      
  contopts = contset();
      
  input.point_on_limitcycle      = [];
  input.odefile                  = [];
  input.ode_parameters           = [];
  input.active_parameter_index   = [];
  input.lower_bound_period       = [];
  input.upper_bound_period       = [];
  input.nMeshIntervals           = [];
  input.subspace_size            = []; 
  % todo: set a default for subspace size, for instance nphases / 2
  input.integration_method  = @ode15s;
  input.integration_options = odeset( ...
    'AbsTol', contopts.time_integration_abs_tol, ...
    'RelTol', contopts.time_integration_rel_tol);
  input.poincare_tolerance       = 1e-2;
  input.show_plot                = false;
  
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
  if input.nMeshIntervals <= 1
    error('nMeshIntervals must be greater than one');
  end
  
  
  initial_continuation_data = init_multiple_shooting_internal(input);
end



