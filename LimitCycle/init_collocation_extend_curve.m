%% initialize a cycle continuation using collocation from a previous run
% The arguments to this function must be specified as name value pairs, for
% instance: (see also Tutorial/fusion/extend_fusion_cycles.m):
%
% load(point_file, 'point');
%
% initial_continuation_data = init_collocation_extend_curve( ...
%   'continuation_state',           point.x, ...
%   'continuation_tangent',         point.v, ... 
%   'odefile',                      odefile,  ...
%   'ode_parameters',               point.parametervalues, ...
%   'active_parameter_index',       3, ...
%   'time_mesh',                    point.timemesh, ...
%   'current_nMeshIntervals',       point.ntst, ...
%   'current_nCollocationPoints',   point.ncol ...
% );
%
%
%% +++++ required arguments +++++
%
%% continuation_state
% The continuation state vector of point from which the continuation should be
% extended. The continuation state vector for continuation of cycles using
% collocation in cl_matcontL contains (current_nMeshIntervals *
% current_nCollocationPoints + 1) * nphase + 2 elements. The first
% (current_nMeshIntervals * current_nCollocationPoints + 1) * nphase points
% contain the coordinates of the points on the cycle. The first nphase elements
% contain the coordinates of the first point, and so on. The first and last
% points are the same. The butlast element of the contiuation state vector
% contains the period, and the last element contains the value of the active
% parameter.
%
%% odefile
% A function handle of the odefile that specifies the system of ODEs. A jacobian
% of the system of ODEs does not need to be specified. If the Jacobian is not
% specified, it will be computed using finite differences ( in the file
% Continuer/cjac.m ). Note that using finite differences in very large systems
% (more than one hundred equations ), might slow down the continuation
% significantly.
%
%% ode_parameters
% the parameter values of the system of ODEs, at which the extension of the
% continuation starts.
%
%% active_parameter_index
% the 1-based index of the parameter that is to be varied during the
% continuation
%
%% time_mesh
% The time_mesh of the cycle from which the extension is started. When using
% collocation, the period of the cycle is subdivided into nMeshIntervals
% intervals. The width of these intervals is adapted, to optimize the
% performance of the continuation. Thus, the time_mesh is variable, and must be
% provided to make the extension of a continuation possible. time_mesh contains
% nMeshIntervals + 1 points. Each interval in further subdivided into
% nCollocationPoints parts to obtain a fine mesh. Within each mesh interval the
% fine mesh is equidistant, and can be obtained using
% Limitcycle/get_fine_mesh.m.
%
%% current_nMeshIntervals
% The number of mesh intervals used in the cycle continuation that is to be
% extended.
%
%% current_nCollocationPoints
% The number of collocation points used in the cycle continuation that is to be
% extended. This number is at least 2 and at most 7.
%
%% +++++ optional arguments +++++
%
%% new_nMeshIntervals
% The number of mesh intervals that is to be used in the extension.
%
%% new_nCollocationPoints
% The number of collocation points that is to be used in the extension.
%
%% continuation_tangent_vector
% The continuation tangent vector at the point from which the continuation
% should be continued. If it is supplied the direction of the continuation can
% umambiguously specified. (by setting the continuer option set_direction to
% false ( see contset.m or the code in contL.m)). If it is not supplied, the
% tangent vector is computed, and the continuation direction is in the direction
% that correspond to an increase in the active_parameter ( if the component of
% the tangent vector is not negligibly small, and if contopts.Backward == false
% ( see contset.m ))


function [continuation_state, continuation_tangent] = ...
                                         init_collocation_extend_curve(varargin)
	% todo: check if nCollocationPoints arugments are between 2 and 7 (see
	% nc_weight.m)
  must_be_specified_by_user = [];
  input.continuation_state          = must_be_specified_by_user;
  input.continuation_tangent        = 0;
  input.odefile                     = must_be_specified_by_user;
  input.ode_parameters              = must_be_specified_by_user;
  input.active_parameter_index      = must_be_specified_by_user;
  input.time_mesh                   = must_be_specified_by_user;
  input.current_nMeshIntervals      = must_be_specified_by_user;
  input.current_nCollocationPoints  = must_be_specified_by_user;
  input.new_nMeshIntervals          = 0; % default value: current_nMeshIntervals
  input.new_nCollocationPoints      = 0; % default value: current_nCollocationPoints
  
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
  
  if input.continuation_tangent == 0
    input.continuation_tangent = [];
  end
  
  if input.new_nMeshIntervals == 0
    input.new_nMeshIntervals = input.current_nMeshIntervals;
  end
  
  if input.new_nCollocationPoints == 0
    input.new_nCollocationPoints = input.current_nCollocationPoints;
  end
  
  [continuation_state, continuation_tangent] = ...
          do_init_collocation_extend_curve(input);
end

function [continuation_state, continuation_tangent] = ...
                                    do_init_collocation_extend_curve(in)

  global lds cds
  
  cds = [];
  lds = [];
  
  n_par = length(in.active_parameter_index);
  if n_par ~= 1 && n_par ~= 2
    error(['One active parameter and the period or 2 active parameters ' ...
           'are needed for limt cycle continuation']);
  end

  odefile_handles = feval(in.odefile);
 
  cds.options = contset();
  
  [max_order, max_order_params] = ...
    find_maximum_order_of_symbolic_derivatives(in.odefile);
  
  cds.options       = contset(cds.options, 'SymDerivative',  max_order);
  cds.options       = contset(cds.options, 'SymDerivativeP', max_order_params);
  cds.symjac        = 1;
  cds.symhess       = 0;
  cds.probfile      = in.odefile;
  cds.P0            = in.ode_parameters;
  cds.ncoo          = length(in.continuation_state) - 1;
  cds.nap           = 1;
  cds.ActiveParams  = in.active_parameter_index;
  cds.usernorm      = odefile_handles{10};
  
  if length(odefile_handles) > 10      
    cds.userfunc = odefile_handles(10:end);
  else
    cds.userfunc = [];
  end
  
  initialize_lds_fields();
  loadOdeFile_into_lds(in.odefile);
  
  nPoint_coordinates = size(in.continuation_state, 1);
  nPoints_on_cycle   = ...
    in.current_nMeshIntervals * in.current_nCollocationPoints + 1;
  
  lds.nphase = round(nPoint_coordinates / nPoints_on_cycle);
  lds.ActiveParams = in.active_parameter_index; % duplicate of cds.Activeparams, todo: carefully remove
  lds.P0           = in.ode_parameters;         % duplicate of cds.Activeparams, todo: carefully remove
  set_ntst_ncol(in.current_nMeshIntervals, ...
                in.current_nCollocationPoints, ...
                in.time_mesh);
  lds.T = in.continuation_state(end-1); % used for fixed period limitcycles
  
  % generate a new mesh and interpolate
  [continuation_state, continuation_tangent] = new_mesh(...
                        in.continuation_state, ...
                        in.continuation_tangent, ...
                        in.new_nMeshIntervals, ...
                        in.new_nCollocationPoints);
  % note: set_ntst_ncol is called from new_mesh
  
  lds.ups = []; 
  lds.vps = [];
  lds.BranchParams=[]; 
end
