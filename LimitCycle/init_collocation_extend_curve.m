function [continuation_state, continuation_tangent] = ...
                                         init_collocation_extend_curve(varargin)
	% todo: check if nCollocationPoints arugments are between 2 and 7 (see
	% nc_weight.m)
  must_be_specified_by_user = [];
  input.continuation_state          = must_be_specified_by_user;
  input.continuation_tangent        = must_be_specified_by_user;
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
