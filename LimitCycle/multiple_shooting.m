% Curve file of cycle continuation of cycles with multiple shooting
%
% The way a cycle is represented as a vector x in cl_matcontL when using
% multiple shooting is as follows. The first cds.n_phases coordinates of x are
% coordinates of a point on the cycle. The next cds.nphases coordinates of x are
% coordinates of the next point on the cycle. The butlast coordinate of x is the
% period of the cycle, and the last coordinate of x is the value of the active
% parameter.
%
% The time mesh is adjusted every few steps. The interval between adjustments
% can be configured using the Adapt settings ( see contset.m ). 
function out = multiple_shooting
  out{1}  = @curve_func;
  out{2}  = @default_processor;
  out{3}  = @options;
  out{4}  = @jacobian;
  out{5}  = [];%@hessians;
  out{6}  = @testfunctions;
  out{7}  = [];%@userf;
  out{8}  = @shooting_process_singularity;
  out{9}  = @cycle_singularity_matrix;
  out{10} = @locate;
  out{11} = @init;
  out{12} = [];%@done;
  out{13} = @adapt;
  out{14} = @curve_CIS_first_point;
  out{15} = @curve_CIS_step;  
end
%-------------------------------------------------------------------------------
% Computes the curve function. If every component of func is zero,  then y_0
% corresponds to a sequence of points on a limit cycle, the "period" is its
% period, and "parameters" are the parameter values of the parameters of the
% system of ODEs for which the cycle exists.
function func = curve_func(varargin)
  global cds
  [y_0, period, parameters] = getComponents(varargin{1});
  y_end = zeros(cds.n_phases * cds.n_mesh_intervals, 1);
  for i = 0 : cds.n_mesh_intervals - 1
    indices = (1:cds.n_phases) + i * cds.n_phases;
    delta_t = period * (cds.mesh(i+2) - cds.mesh(i+1));
    y_end(indices) = NP_shoot(y_0(indices), delta_t, parameters);
  end
  r = zeros(cds.n_phases * cds.n_mesh_intervals,1); % residuals
  for i = 0 : cds.n_mesh_intervals - 2
    indices1 = (1:cds.n_phases) +  i      * cds.n_phases;
    indices2 = (1:cds.n_phases) + (i + 1) * cds.n_phases;
    r(indices1) = y_end(indices1) - y_0(indices2);
  end
  r((1:cds.n_phases) + (cds.n_mesh_intervals-1) * cds.n_phases) = ... 
          y_end(end-cds.n_phases+1:end) - y_0(1:cds.n_phases);
  func = [r;(y_0(1:cds.n_phases) - cds.previous_phases)' * cds.previous_dydt_0]; 
end
%-------------------------------------------------------------------------------
% Computes the Jacobian matrix of the curvefunction at evaluated at varargin{1}.
function jacobian = jacobian(varargin)
  global cds
  M = cds.n_mesh_intervals * cds.n_phases;
  n_phases = cds.n_phases;
  m = cds.n_mesh_intervals;
  % compute number of nonzero's (nnz)
  nnz = m * n_phases^2 + m * n_phases; % main blocks
  nnz = nnz + 2 * (M + 1);           % last 2 columns
  nnz = nnz + M + 2;                 % bottom row
  jacobian = spalloc(M + 1, M + 2,nnz);
  [y_0, period, parameters] = getComponents(varargin{1});
  y_end = zeros(M,1);
  for i = 0 : cds.n_mesh_intervals - 1
    indices = (1:cds.n_phases) + i * cds.n_phases;
    [y_end(indices), jacobian(indices,indices)] = ...
            compute_monodromy(y_0(indices), i + 1, period, parameters);
  end
  
  for i = 0 : cds.n_mesh_intervals - 2
    indices1 = (1:cds.n_phases) + i     * cds.n_phases;
    indices2 = (1:cds.n_phases) + (i+1) * cds.n_phases;
    jacobian(indices1, indices2) = - eye(cds.n_phases); %#ok<*SPRIX>
    % This jacobian function is used for testing only, therefore, we ignore the
    % "this sparse indexing operation is likely to be slow"-warning.
  end
  jacobian((1:cds.n_phases) + (cds.n_mesh_intervals-1) * cds.n_phases, ...
            1:cds.n_phases) = - eye(cds.n_phases);
  for i = 0 : cds.n_mesh_intervals - 1
    % compute d_y_d_T
    indices = (1:cds.n_phases) + i * cds.n_phases;
    jacobian(indices,M+1) = cds.dydt_ode(0, y_end(indices), parameters{:}) * ...
      (cds.mesh(i+2) - cds.mesh(i+1));
  end
  for i = 0 : cds.n_mesh_intervals - 1
    % compute d_y_d_p
    indices = (1:cds.n_phases) + i * cds.n_phases;
    delta_t = period * (cds.mesh(i+2) - cds.mesh(i+1));
    jacobian(indices,M+2) = ...
       NP_compute_d_phi_d_p(y_0(indices), delta_t, parameters);
  end
  % specify d_s_d_y
  jacobian(M+1, 1:cds.n_phases) = cds.previous_dydt_0';
  % specify d_s_d_T = jacobian(N+1,N+1) = 0;
  % specify d_s_d_p = jacobian(N+1,N+2) = 0;
end
%-------------------------------------------------------------------------------
% computes d_phi over d_x where phi is the solution of the initial value problem
% phi' = f(phi),   phi(0) = x
% the result is a square Jacobian matrix of size cds.n_phases
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
  n_phases = cds.n_phases;
  f = @(t, y) dydt_monodromy_full(t, y, parameters);
  integration_opt = odeset(...
    'AbsTol',      contopts.integration_abs_tol,    ...
    'RelTol',      contopts.integration_rel_tol     ... % todo add JPattern   
  );
  time_interval = period * [cds.mesh(mesh_index)  cds.mesh(mesh_index+1)];
  x_with_monodromy = [x; reshape(eye(n_phases),[n_phases^2 1])];
  [~, trajectory] = cds.integrator( ...
                        f, time_interval, x_with_monodromy, integration_opt);
  y_end = trajectory(end,1:n_phases)';
  monodromy = trajectory(end,n_phases+1:end);
  monodromy = reshape(monodromy, [n_phases n_phases]);
end
%-------------------------------------------------------------------------------
function dydt_combined = dydt_monodromy_full(t, y, parameters)
  global cds
  y_ode          = y(1:cds.n_phases);
  jacobian       = cds.jacobian_ode(t, y_ode, parameters{:});
  monodromy      = reshape(y(cds.n_phases+1:end), cds.n_phases, cds.n_phases);
  dydt_monodromy = jacobian * monodromy;
  dydt_combined  = [
      cds.dydt_ode(t, y_ode, parameters{:}); 
      dydt_monodromy(:)
  ];
end
%------------------------------------------------------------------------------- 
function [y_end, monodromy] = ...
                   monodromy_column_by_column(x, mesh_index, period, parameters)
  global cds contopts;
  % compute trajectory of the cycle
  integration_opt_1 = odeset(...
    'AbsTol',       contopts.integration_abs_tol, ...
    'RelTol',       contopts.integration_rel_tol, ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode, t, y, parameters{:}));
  f                 = @(t, y) cds.dydt_ode(t, y, parameters{:});
  time_interval     = period * [cds.mesh(mesh_index) cds.mesh(mesh_index+1)];
  cycle             = cds.integrator(f, time_interval, x, integration_opt_1);
  y_end             = deval(cycle, period);
  
  % compute the monodromy matrix
  monodromy         = eye(cds.n_phases);
  integration_opt_2 = odeset( ...
    'AbsTol',       contopts.integration_abs_tol, ...
    'RelTol',       contopts.integration_rel_tol, ...
    'Jacobian',     @(t,y) cds.jacobian_ode(t, deval(cycle, t), parameters{:}));
  f = @(t, y) cds.jacobian_ode(t, deval(cycle,t), parameters{:}) * y;
  
  for i = 1 : cds.n_phases
    [~, monodromy_map_trajectory] = cds.integrator(...
            f, time_interval, monodromy(:,i), integration_opt_2);
          
    monodromy(:,i) = monodromy_map_trajectory(end,:);
  end 
end
%-------------------------------------------------------------------------------
% Test functions are used for detecting AND locating singularities by bisection.
% When detecting ids_testf_requested will be cds.ActTest, and when locating
% ids_testf_requested will contain only those ids of the testfunctions relevant
% to the bifurcation that is being located.
function [out, failed] = testfunctions(ids_testf_requested, x0, v, ~) 
  global cds
  
  failed = false;
  
  const = Constants;
  
  if any(ismember([const.BPC_id const.PD_id const.NS_id], ids_testf_requested))
    update_multipliers_if_needed(x0)
  end
  
  out = cycle_testfunctions(ids_testf_requested, cds.multipliers, v);
end
%-------------------------------------------------------------------------------
function init(~,~); end
%-------------------------------------------------------------------------------
function point = default_processor(varargin)
  global cds
  point = varargin{1};
  
  point = adjust_basis_size(point);

  [y, ~, parameter_values] = getComponents(point.x);
  cds.previous_phases      = y(1:cds.n_phases);
  cds.previous_dydt_0      = cds.dydt_ode(0, cds.previous_phases, ...
                                                           parameter_values{:});

 
  update_multipliers_if_needed(point.x);
  
  point.multipliers      = cds.multipliers;
  point.parameter_values = cell2mat(parameter_values);
  point.n_mesh_intervals = cds.n_mesh_intervals;
  point.time_mesh        = cds.mesh;
  point.n_phases         = cds.n_phases;
  savePoint(point, varargin{2:end});
end
%-------------------------------------------------------------------------------
function options
end
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
  indices = (1:cds.n_phases) + cds.n_phases;
  for i = 2 : cds.n_mesh_intervals
    x(indices) = interp1( ...
            cds.t_cycle, cds.y_cycle, period * cds.mesh(i), 'spline');
    indices = indices + cds.n_phases;
  end
  point.R = max(abs(curve_func(x)));
  print_diag(4, 'new_time_mesh:');
  print_diag(4, ' %.4f', cds.mesh );
  print_diag(4, '\n');
  print_diag(4, 'curve_function new time mesh: %.3e\n', point.R);
end
%-------------------------------------------------------------------------------
function point = adjust_basis_size(point)
  global cds contopts;
  basis_size_changed = false;
  
  if contopts.NewtonPicard
    update_multipliers_if_needed(point.x)
    if abs(cds.multipliers(end)) > contopts.basis_grow_threshold
      basis_size_changed = true;
      print_diag(2, 'expanding basis\n');
      nMults_to_compute = cds.preferred_basis_size + 10;
      nMults_to_compute = min(nMults_to_compute, cds.n_phases);
      cds.multipliersX = point.x;
      cds.multipliers = NP_MS_compute_multipliers(...
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
    test_function_labels = {'BPC', 'PD', 'LPC', 'NS'};
    print_diag(1,'testfunctions: [')
    for i=1:length(point.tvals)
      print_diag(1,' %s:',test_function_labels{i});
      print_diag(1,'%+.5e',point.tvals(i))
    end
    print_diag(1,']\n')
  end
end
%-------------------------------------------------------------------------------
function [y,period,parameters] = getComponents(x)
  global cds
  y                            = x(1 : cds.n_phases * cds.n_mesh_intervals);
  period                       = x(end-1);
  parameter_value              = x(end);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = parameter_value;
  parameters                   = num2cell(parameters);
end
%-------------------------------------------------------------------------------
function mesh_points = find_mesh_points_multiple_shooting(x, fineness)
  if nargin == 1
    fineness = 10000;
  end
  global cds contopts
  int_opt = odeset( ...
    'AbsTol', contopts.integration_abs_tol, ...
    'RelTol', contopts.integration_rel_tol  ...
  );
  [~, period, parameters] = getComponents(x);
  NP_MS_compute_stiched_orbit(x);
  cycle_gradient_norm = @(t,x) norm(cds.dydt_ode(...
      t, interp1(cds.t_cycle,cds.y_cycle,t,'spline'), parameters{:}));
  
  time_points = linspace(0, period, fineness * cds.n_mesh_intervals);
  [t,x] = ode15s(cycle_gradient_norm, time_points, 0, int_opt);
  cycle_gradient_integral = x(end);
  x = mod(x, cycle_gradient_integral / cds.n_mesh_intervals);
  mesh_points = zeros(cds.n_mesh_intervals + 1, 1);
  mesh_points_index = 2;
  for i=2:size(x)
    if x(i-1) > x(i)
      mesh_points(mesh_points_index) = t(i);
      mesh_points_index = mesh_points_index + 1;
    end
  end
  mesh_points = mesh_points / period;
  % if the gradient increases more than cycle_gradient_integral /
  % cds.n_mesh_intervals per time interval ( i.e. consecutive values of t, then
  % the the for loop might not produce enough mesh points. Hence, in this case
  % we must reduce the time interval and try again.
  if (mesh_points(end) == 0)
    mesh_points = find_mesh_points_multiple_shooting(x, fineness * 10);
  end
end
%-------------------------------------------------------------------------------
function update_multipliers_if_needed(x)
  global cds
  if ~ isfield(cds,'multipliersX') || all(cds.multipliersX ~= x)
    cds.multipliersX = x;
    cds.multipliers = NP_MS_compute_multipliers(x, cds.preferred_basis_size);
  end
end
%-------------------------------------------------------------------------------
function p_out = locate(id, p1, p2) %#ok<INUSD,STOUT>
  switch id   
    otherwise
      error('No locator defined for singularity %d', id);
  end
end
%-------------------------------------------------------------------------------

