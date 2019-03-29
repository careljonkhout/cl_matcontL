function out = multiple_shooting
%
% Curve file of cycle continuation with multiple shooting
%
    out{1}  = @curve_func;
    out{2}  = @defaultprocessor;
    out{3}  = @options;
    out{4}  = @jacobian;
    out{5}  = [];%@hessians;
    out{6}  = [];%@testf;
    out{7}  = [];%@userf;
    out{8}  = @process;
    out{9}  = [];%@singmat;
    out{10} = [];%@locate;
    out{11} = @init;
    out{12} = [];%@done;
    out{13} = @adapt;
    out{14} = @curve_CIS_first_point;
    out{15} = @curve_CIS_step;
%---------------------------------------------------------  
end

function func = curve_func(varargin)
  global cds
  [y_0, period, parameters] = getComponents(varargin{1});
  y_end = zeros(cds.nphases * cds.nMeshPoints,1);
  for i=0:cds.nMeshPoints-1
    indices = (1:cds.nphases) + i * cds.nphases;
    y_end(indices) = ...
      shoot(y_0(indices), i+1, period, parameters);
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
%---------------------
  
function x_end = shoot(x, mesh_index, period, parameters)
  global cds contopts
  f =@(t, y) cds.dydt_ode(t, y, parameters{:});
  integration_opt = odeset(...
    'AbsTol',       contopts.integration_abs_tol,    ...
    'RelTol',       contopts.integration_rel_tol,    ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  time_interval = period * [cds.mesh(mesh_index) cds.mesh(mesh_index+1)];
  [~, trajectory] = ode15s(f, time_interval, x, integration_opt);
  x_end = trajectory(end,:)';
end

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
    jacobian(indices,M+2) = compute_d_phi_d_p(...
      y_0(indices), i+1, period, parameters);
  end
  % specify d_s_d_y
  jacobian(M+1,1:cds.nphases) = cds.previous_dydt_0';
  % specify d_s_d_T = jacobian(N+1,N+1) = 0;
  % specify d_s_d_p = jacobian(N+1,N+2) = 0;
end
  
function d_phi_d_p = compute_d_phi_d_p(x, mesh_index, period, parameters)
  global cds
  ap = cds.ActiveParams;
  h = 1e-6;
  parameters{ap} = parameters{ap} - h;
  phi_1 = shoot(x, mesh_index, period, parameters);
  parameters{ap} = parameters{ap} + 2*h;
  phi_2 = shoot(x, mesh_index, period, parameters);
  d_phi_d_p = (phi_2 - phi_1)/h/2;  
end

function [y_end, monodromy] = compute_monodromy(x, mesh_index, period, parameters)
  global cds
  if ~ cds.options.PartitionMonodromy
    [y_end, monodromy] = monodromy_full(x, mesh_index, period, parameters);
  else
    [y_end, monodromy] = ...
      monodromy_column_by_column(x, mesh_index, period, parameters);
  end
end
    
    
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

function init(~,~)
end

%----------------------------------------------------------
function point = defaultprocessor(varargin)
  global cds contopts
  point = varargin{1};
  point.mesh = cds.mesh;
  if contopts.Multipliers
    point.multipliers = ...
      NewtonPicard.MultipleShooting.compute_multipliers(point.x);
  end
  [y, ~, parameters]      = getComponents(point.x);
  cds.previous_phases     = y(1:cds.nphases);
  cds.previous_dydt_0     = cds.dydt_ode(0, cds.previous_phases, parameters{:});

  savePoint(point);
end

function options
end

function CISdata = curve_CIS_first_point(x) %#ok<INUSD>
  CISdata = 1;
end

function CISdata = curve_CIS_step(x, CISdata_in) %#ok<INUSD>
  CISdata = 1;
end

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

function [y,period,parameters] = getComponents(x)
  global cds
  y                            = x(1:cds.nphases*cds.nMeshPoints);
  period                       = x(end-1);
  parameter_value              = x(end);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = parameter_value;
  parameters                   = num2cell(parameters);
end

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

function [failed, s] = process(si, singpoint, s) %#ok<INUSL>
  % to be implemented later
  failed = false;
end

