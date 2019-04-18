% Curve file of cycle continuation with multiple shooting
function out = multiple_shooting
  out{1}  = @curve_func;
  out{2}  = @defaultprocessor;
  out{3}  = @options;
  out{4}  = @jacobian;
  out{5}  = [];%@hessians;
  out{6}  = @testfunctions;
  out{7}  = [];%@userf;
  out{8}  = @process_singularity;
  out{9}  = @singularity_matrix;
  out{10} = [];%@locate;
  out{11} = @init;
  out{12} = [];%@done;
  out{13} = @adapt;
  out{14} = @curve_CIS_first_point;
  out{15} = @curve_CIS_step;  
end
%-------------------------------------------------------------------------------
function func = curve_func(varargin)
  global cds
  [y_0, period, parameters] = getComponents(varargin{1});
  y_end = zeros(cds.nphases * cds.nMeshPoints,1);
  for i=0:cds.nMeshPoints-1
    indices = (1:cds.nphases) + i * cds.nphases;
    delta_t = period * [cds.mesh(i+1) cds.mesh(i+2)];
    y_end(indices) = NewtonPicard.shoot(y_0(indices), delta_t, parameters);
  end
  r = zeros(cds.nphases * cds.nMeshPoints,1); % residuals
  for i=0:cds.nMeshPoints-2
    indices1 = (1:cds.nphases) + i     * cds.nphases;
    indices2 = (1:cds.nphases) + (i+1) * cds.nphases;
    r(indices1) = y_end(indices1) - y_0(indices2);
  end
  r((1:cds.nphases)+(cds.nMeshPoints-1)*cds.nphases) ... 
    = y_end(end-cds.nphases+1:end) - y_0(1:cds.nphases);
%  cds.cycle_trajectory = ode15s(...
%    @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
%    linspace(0, period, cds.nDiscretizationPoints), ...
%    phases_0, integration_opt);
  func = [r
          (y_0(1:cds.nphases) - cds.previous_phases)' * cds.previous_dydt_0 ]; 
end
%-------------------------------------------------------------------------------
function jacobian = jacobian(varargin)
  global cds
  M = cds.nMeshPoints * cds.nphases;
  nphases = cds.nphases;
  m = cds.nMeshPoints;
  nnz = m * nphases^2 + m * nphases; % main blocks
  nnz = nnz + 2 * (M + 1);           % last 2 columns
  nnz = nnz + M + 2;                 % bottom row
  jacobian = spalloc(M + 1, M + 2,nnz);
  [y_0, period, parameters] = getComponents(varargin{1});
  y_end = zeros(M,1);
  for i=0:cds.nMeshPoints-1
    indices = (1:cds.nphases) + i * cds.nphases;
    [y_end(indices), jacobian(indices,indices)] = ...
      compute_monodromy(y_0(indices), i+1, period, parameters);
  end
  
  for i=0:cds.nMeshPoints-2
    indices1 = (1:cds.nphases) + i     * cds.nphases;
    indices2 = (1:cds.nphases) + (i+1) * cds.nphases;
    jacobian(indices1, indices2) = - eye(cds.nphases); %#ok<*SPRIX>
    % This jacobian function is used for testing only, therefore, we ignore the
    % "this sparse indexing operation is likely to be slow"-warning.
  end
  jacobian((1:cds.nphases) + (cds.nMeshPoints-1) * cds.nphases, ...
            1:cds.nphases) = - eye(cds.nphases);
  for i=0:cds.nMeshPoints-1
    % compute d_y_d_T
    indices = (1:cds.nphases) + i * cds.nphases;
    jacobian(indices,M+1) = cds.dydt_ode(0, y_end(indices), parameters{:}) * ...
      (cds.mesh(i+2) - cds.mesh(i+1));
  end
  for i=0:cds.nMeshPoints-1
    % compute d_y_d_p
    indices = (1:cds.nphases) + i * cds.nphases;
    delta_t = period * [cds.mesh(i+1) cds.mesh(i+2)];
    jacobian(indices,M+2) = ...
       NewtonPicard.compute_d_phi_d_p(y_0(indices), delta_t, parameters);
  end
  % specify d_s_d_y
  jacobian(M+1,1:cds.nphases) = cds.previous_dydt_0';
  % specify d_s_d_T = jacobian(N+1,N+1) = 0;
  % specify d_s_d_p = jacobian(N+1,N+2) = 0;
end
%-------------------------------------------------------------------------------
function [y_end, monodromy] = ...
                            compute_monodromy(x, mesh_index, period, parameters)
  global cds
  if ~ cds.options.PartitionMonodromy
    [y_end, monodromy] = monodromy_full(x, mesh_index, period, parameters);
  else
    [y_end, monodromy] = ...
      monodromy_column_by_column(x, mesh_index, period, parameters);
  end
end
%-------------------------------------------------------------------------------
function [y_end, monodromy] = monodromy_full(x, mesh_index, period, parameters)
  global cds contopts
  nphases = cds.nphases;
  f =@(t, y) dydt_monodromy_full(t, y, parameters);
  integration_opt = odeset(...
    'AbsTol',      contopts.integration_abs_tol,    ...
    'RelTol',      contopts.integration_rel_tol     ... % todo add JPattern   
  );
  time_interval = period * [cds.mesh(mesh_index)  cds.mesh(mesh_index+1)];
  x_with_monodromy = [x; reshape(eye(nphases),[nphases^2 1])];
  [~, trajectory] = ode15s(f, time_interval, x_with_monodromy, integration_opt);
  y_end = trajectory(end,1:nphases)';
  monodromy = trajectory(end,nphases+1:end);
  monodromy = reshape(monodromy, [nphases nphases]);
end
%-------------------------------------------------------------------------------
function dydt_mon = dydt_monodromy_full(t, y, parameters)
  global cds
  y_ode = y(1:cds.nphases);
  
  y_mon = reshape(y(cds.nphases+1:end),cds.nphases,cds.nphases);
  dydt_mon = [
      cds.dydt_ode(t, y_ode, parameters{:}); 
      reshape( ...
        cds.jacobian_ode(t, y_ode, parameters{:}) * y_mon, ...
        [cds.nphases^2 1]) 
  ];
end
%------------------------------------------------------------------------------- 
function [y_end, monodromy] = ...
                   monodromy_column_by_column(x, mesh_index, period, parameters)
  global cds contopts;
  integration_opt = odeset(...
    'AbsTol',       contopts.integration_abs_tol,    ...
    'RelTol',       contopts.integration_rel_tol,    ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  f = @(t, y) cds.dydt_ode(t, y, parameters{:});
  time_interval = period * [cds.mesh(mesh_index) cds.mesh(mesh_index+1)];
  cycle         = ode15s(f, time_interval, x, integration_opt);
  y_end         = deval(cycle,period);
  monodromy     = eye(cds.nphases);
  integration_opt = odeset(integration_opt, 'Jacobian', ...
    @(t,y) feval(cds.jacobian_ode, t, deval(cycle,t), parameters{:}));
  f = @(t, y) cds.jacobian_ode(t, deval(cycle,t), parameters{:}) * y;
  
  parfor i=1:cds.nphases
    print_diag(1,'computing column %d of monodromy matrix ',i);
    [~, monodromy_map_trajectory] = ode15s(...
      f, time_interval, monodromy(:,i), integration_opt);
    monodromy(:,i) = monodromy_map_trajectory(end,:);
  end 
end
%-------------------------------------------------------------------------------
% Test functions are used for detecting AND location singularities by bisection.
% When detecting ids_testf_requested will be cds.ActTest, and when locating
% ids_testf_requested will contain only those ids of the testfunctions relevant
% to the bifurcation that is being located.
function [out, failed] = testfunctions(ids_testf_requested, x0, v, ~) 
  % unused arguments are v and CISdata
  global cds contopts
  
  failed = false;
  
  const = Constants;
  
  if any(ismember([const.BPC_id const.PD_id const.NS_id], ids_testf_requested))
    update_multipliers_if_needed(x0)
  end
  if any(ismember(const.BPC_id, ids_testf_requested))
    mults = cds.multipliers;
    distance_to_one = abs(mults - 1);
    [~,index]       = min(distance_to_one);
    mults(index) = 0;
    out(const.BPC_id) = real(prod(mults - ones(size(mults))));
  end
  if ismember(const.PD_id, ids_testf_requested)
    % real is needed to ensure that small complex parts induced by roundoff
    % errors do not make the result complex.
    % A complex valued test function would cause false positives when detecting 
    % bifurcations.
    out(const.PD_id) = real(prod(cds.multipliers + ones(size(cds.multipliers))));
  end
  if ismember(const.LPC_id, ids_testf_requested)
    out(const.LPC_id) = v(end);
  end
  if any(ismember(const.NS_id, ids_testf_requested))
    mults                  = cds.multipliers;
    threshold              = contopts.real_v_complex_threshold;
    complex_mults          = mults(abs(imag(mults)) > threshold);
    unstable_complex_mults = complex_mults(abs(complex_mults) > 1); 
    out(const.NS_id)       = length(unstable_complex_mults);
  end
  
end
%-------------------------------------------------------------------------------
% defines which changes in testfunctions correspond to which singularity type
% 0 == require sign-change
% 1 == require sign-non-change
% 2 == require change
% anything else: no requirement
% columns correspond to testfunctions
% rows correspond to singularities BPC, PD, LPC, and NS respectively
function [S,L] = singularity_matrix

  
  sign_change    = Constants.sign_change;
  sign_constant  = Constants.sign_constant;
  value_change   = Constants.value_change;
  o              = Constants.ignore; 

  S = [ 
  % BPC_testfunc PD_testfunc  LPC_testfunc   NS_testfunc
    sign_change  o            sign_constant  o            % BPC    
    o            sign_change  o              o            % PD   
    o            o            sign_change    o            % LPC   
    o            o            o              value_change % NS
  ];
  L = [ 'BPC';'PD '; 'LPC'; 'NS ' ];
end
%-------------------------------------------------------------------------------
% After a singularity is detected and located, the contL calls this function.
% Usually a normal form coefficient (nfc) is computed, via a function call from
% this function. However, for continuation of cycles by single shooting or
% multiple shooting, the computation of normal form coefficients is not
% implemented yet.
function [failed,s] = process_singularity(id,point,s)
  global cds
  x = point.x;
  switch id
  case 1
    % note: as of march 2019 detection for BPC not yet implemented 
    format_string = 'Branch Point cycle(period = %e, parameter = %e)\n'; 
    print_diag(0, format_string, x(end-1), x(end));
    s.msg  = sprintf('Branch Point cycle'); 
  case 2
    format_string = 'Period Doubling (period = %e, parameter = %e)\n';
    print_diag(0, format_string, x(end-1), x(end));
    s.msg  = 'Period Doubling';
  case 3
    s.msg = 'Limit point cycle';
    format_string = 'Limit point cycle (period = %e, parameter = %e)\n';
    print_diag(0, format_string, x(end-1), x(end));
  case 4
    d = cds.multipliers;
    smallest_sum = Inf;
    for jk=1:cds.nphases-1
      [val,idx] = min(abs(d(jk+1:cds.nphases)*d(jk)-1));
      if val < smallest_sum
        idx2 = jk+idx;
        smallest_sum = val;
      end
    end
    singularity_is_neutral_saddle = imag(d(idx2)) == 0;
    if singularity_is_neutral_saddle
      s.msg = 'Neutral saddle cycle';
      format_string = 'Neutral Saddle Cycle (period = %e, parameter = %e)\n';
      % A neutral saddle is not really a bifurcation, therefore we use priority 
      % 1 instead of 0, so that it is only logged if 
      % contopts.contL_DiagnosticsLevel is set higher than the default value
      % which is zero.
      print_diag(1, format_string, x(end-1), x(end));
    else
      s.msg = 'Neimark Sacker';
      format_string = 'Neimark-Sacker (period = %e, parameter = %e)\n';
      print_diag(0, format_string, x(end-1) ,x(end));
      %print_diag(0, 'Normal form coefficient = %d\n', s.data.nscoefficient);
    end
  end
  failed = 0;
end  
%-------------------------------------------------------------------------------
function init(~,~); end
%-------------------------------------------------------------------------------
function point = defaultprocessor(varargin)
  global cds
  point = varargin{1};
  point.mesh = cds.mesh;

  [y, ~, parameters]      = getComponents(point.x);
  cds.previous_phases     = y(1:cds.nphases);
  cds.previous_dydt_0     = cds.dydt_ode(0, cds.previous_phases, parameters{:});

  adjust_basis_size(point);
  savePoint(point);
end
%-------------------------------------------------------------------------------
function options; end
%-------------------------------------------------------------------------------
function CISdata = curve_CIS_first_point(~) % unused argument is x 
  CISdata = 1;
end
%-------------------------------------------------------------------------------
function CISdata = curve_CIS_step(~,~) 
  % unused arguments are x and CIS_data_in
  CISdata = 1;
end
%-------------------------------------------------------------------------------
function [has_changed, x, v, CISData] = adapt(varargin)
  global cds
  has_changed = true;
  x       = varargin{1};
  v       = varargin{2};
  CISData = varargin{3};
  period  = x(end-1);
  print_diag(3,'curve_function with old time mesh: %.3e\n', ...
    max(abs(curve_func(x))));
  % We reposition the mesh points in time:
  cds.mesh = find_mesh_points_multiple_shooting(x);
  
  % We update coordinates of the points on the cycle to correspond to the new
  % time mesh. Note that the first point does not change.
  indices = (1:cds.nphases) + cds.nphases;
  for i=2:cds.nMeshPoints
    x(indices) = interp1( ...
      cds.t_cycle, cds.y_cycle,period * cds.mesh(i), 'spline');
    indices = indices + cds.nphases;
  end
  print_diag(4, 'new_time_mesh:');
  print_diag(3, ' %.4f', cds.mesh );
  print_diag(4, '\n');
  print_diag(4, 'curve_function new time mesh: %.3e\n', ...
    max(abs(curve_func(x))));
  %x = NewtonPicard.MultipleShooting.corrections_without_tangent(x);
  %print_diag(3,'curve_function new time mesh: %.3e\n', ...
  %  max(abs(curve_func(x))));
  
end
%-------------------------------------------------------------------------------
function adjust_basis_size(point)
  global cds contopts;
  basis_size_changed = false;
  
  if contopts.NewtonPicard
    update_multipliers_if_needed(point.x)
    if abs(cds.multipliers(end)) > contopts.basis_grow_threshold
      basis_size_changed = true;
      print_diag(2, 'expanding basis\n');
      nMults_to_compute = cds.preferred_basis_size + 10;
      nMults_to_compute = min(nMults_to_compute, cds.nphases);
      cds.multipliersX = point.x;
      cds.multipliers = NewtonPicard.MultipleShooting.compute_multipliers(...
        point.x, nMults_to_compute);
      
      i = length(cds.multipliers);
     
      while abs(cds.multipliers(i)) < contopts.basis_grow_threshold / 3
        i = i - 1;
      end
      if i < length(cds.multipliers)
        i = i + 1;
      end
      cds.preferred_basis_size = i;
      cds.p                    = i;
      cds.multipliers = cds.multipliers(1:i);
    elseif abs(cds.multipliers(end)) < contopts.basis_shrink_threshold
      basis_size_changed = true;
      i = length(cds.multipliers);
      while abs(cds.multipliers(i)) < contopts.basis_shrink_threshold
        i = i - 1;
      end
      cds.preferred_basis_size = i;
      cds.p = i;
      cds.multipliers = cds.multipliers(1:i);
    end
  end
  
  if contopts.Multipliers || contopts.NewtonPicard
    update_multipliers_if_needed(point.x)
    point.multipliers = cds.multipliers;
  end
  
  if basis_size_changed
    point.tvals = testfunctions(cds.ActTest,point.x,point.v,[]);
    print_diag(1,'testfunctions: [')
    print_diag(1,' %+.5e',point.tvals)
    print_diag(1,']\n')
  end
end
%-------------------------------------------------------------------------------
function [y,period,parameters] = getComponents(x)
  global cds
  y                            = x(1:cds.nphases*cds.nMeshPoints);
  period                       = x(end-1);
  parameter_value              = x(end);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = parameter_value;
  parameters                   = num2cell(parameters);
end
%-------------------------------------------------------------------------------
function mesh_points = find_mesh_points_multiple_shooting(x)
  global cds contopts
  int_opt = odeset( ...
    'AbsTol', contopts.integration_abs_tol, ...
    'RelTol', contopts.integration_rel_tol  ...
  );
  [~, period, parameters] = getComponents(x);
  NewtonPicard.MultipleShooting.compute_stiched_orbit(x);
  cycle_gradient_norm = @(t,x) norm(cds.dydt_ode(...
      t, interp1(cds.t_cycle,cds.y_cycle,t,'spline'), parameters{:}));
  
  time_points = linspace(0,period,100*cds.nMeshPoints);
  [t,x] = ode45(cycle_gradient_norm, time_points, 0, int_opt);
  cycle_gradient_integral = x(end);
  x = mod(x, cycle_gradient_integral / cds.nMeshPoints);
  mesh_points = zeros(cds.nMeshPoints + 1, 1);
  mesh_points_index = 2;
  for i=2:size(x)
    if x(i-1) > x(i)
      mesh_points(mesh_points_index) = t(i);
      mesh_points_index = mesh_points_index + 1;
    end
  end
  mesh_points = mesh_points / period;
end
%-------------------------------------------------------------------------------
function update_multipliers_if_needed(x)
  global cds
  if ~ isfield(cds,'multipliersX') || all(cds.multipliersX ~= x)
    cds.multipliersX = x;
    cds.multipliers = NewtonPicard.MultipleShooting.compute_multipliers(x, ...
                       cds.preferred_basis_size);
  end
end
%-------------------------------------------------------------------------------
