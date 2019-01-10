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
    out{8}  = [];%@process;
    out{9}  = [];%@singmat;
    out{10} = [];%@locate;
    out{11} = @init;
    out{12} = [];%@done;
    out{13} = @adapt;
    out{14} = @curve_CIS_first_point;
    out{15} = @curve_CIS_step;
%---------------------------------------------------------  
function func = curve_func(varargin)
  global cds
  [y_0,period,parameters] = getComponents(varargin{1});
  y_end = zeros(cds.nphases * cds.nShootingPoints,1);
  for i=0:cds.nShootingPoints-1
    indices = (1:cds.nphases) + i * cds.nphases;
    y_end(indices) = ...
      shoot(y_0(indices), period / cds.nShootingPoints, parameters);
  end
  r = zeros(cds.nphases * cds.nShootingPoints,1); % residuals
  for i=0:cds.nShootingPoints-2
    indices1 = (1:cds.nphases) + i     * cds.nphases;
    indices2 = (1:cds.nphases) + (i+1) * cds.nphases;
    r(indices1) = y_end(indices1) - y_0(indices2);
  end
  r((1:cds.nphases)+(cds.nShootingPoints-1)*cds.nphases) ... 
    = y_end(end-cds.nphases+1:end) - y_0(1:cds.nphases);
%  cds.cycle_trajectory = ode15s(...
%    @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
%    linspace(0, period, cds.nDiscretizationPoints), ...
%    phases_0, integration_opt);
  func = [r
          (y_0(1:cds.nphases) - cds.previous_phases)' * cds.previous_dydt_0 ]; 
%---------------------
  
function x_end = shoot(x, period, parameters)
  global cds
  f =@(t, y) cds.dydt_ode(t, y, parameters{:});
  integration_opt = odeset(...
    'AbsTol',      1e-10,    ...
    'RelTol',      1e-10,    ...
    'BDF',         'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1,      ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  [~, trajectory] = ode15s(f, [0 period], x, integration_opt);
  x_end = trajectory(end,:)';

function jacobian = jacobian(varargin)
  global cds
  M = cds.nShootingPoints * cds.nphases;
  jacobian = zeros(M + 1, M + 2);
  [y_0,period,parameters] = getComponents(varargin{1});
  y_end = zeros(M,1);
  for i=0:cds.nShootingPoints-1
    indices = (1:cds.nphases) + i * cds.nphases;
    [y_end(indices), jacobian(indices,indices)] = ...
      compute_monodromy(y_0(indices), period / cds.nShootingPoints, parameters);
  end
  
  for i=0:cds.nShootingPoints-2
    indices1 = (1:cds.nphases) + i     * cds.nphases;
    indices2 = (1:cds.nphases) + (i+1) * cds.nphases;
    jacobian(indices1, indices2) = - eye(cds.nphases);
  end
  jacobian((1:cds.nphases) + (cds.nShootingPoints-1) * cds.nphases, ...
            1:cds.nphases) = - eye(cds.nphases);
  for i=0:cds.nShootingPoints-1
    % compute d_y_d_T
    indices = (1:cds.nphases) + i * cds.nphases;
    jacobian(indices,M+1) = cds.dydt_ode(0, y_end(indices), parameters{:});
  end
  for i=0:cds.nShootingPoints-1
    % compute d_y_d_p
    indices = (1:cds.nphases) + i * cds.nphases;
    jacobian(indices,M+2) = compute_d_phi_d_p(...
      y_0(indices), period / cds.nShootingPoints, parameters);
  end
  % specify d_s_d_y
  jacobian(M+1,1:cds.nphases) = cds.previous_dydt_0';
  % specify d_s_d_T = jacobian(N+1,N+1) = 0;
  % specify d_s_d_p = jacobian(N+1,N+2) = 0;
  
  
function d_phi_d_p = compute_d_phi_d_p(x, period, parameters)
  global cds
  ap = cds.ActiveParams;
  h = 1e-6;
  parameters{ap} = parameters{ap} - h;
  phi_1 = shoot(x, period, parameters);
  parameters{ap} = parameters{ap} + 2*h;
  phi_2 = shoot(x, period, parameters);
  d_phi_d_p = (phi_2 - phi_1)/h/2;  
  
function [y_end, monodromy] = compute_monodromy(x, period, parameters)
  global cds
  if ~ cds.options.PartitionMonodromy
    [y_end, monodromy] = monodromy_full(x, period, parameters);
  else
    [y_end, monodromy] = monodromy_column_by_column(x, period, parameters);
  end
    
    
function [y_end, monodromy] = monodromy_full(x, period, parameters)
  global cds
  nphases = cds.nphases;
  f =@(t, y) dydt_monodromy_full(t, y, parameters);
  integration_opt = odeset(...
    'AbsTol',      1e-10,    ...
    'RelTol',      1e-10,    ...
    'BDF',         'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1 ... % todo add jacobian   
  );

  x_with_monodromy = [x; reshape(eye(nphases),[nphases^2 1])];
  [~, trajectory] = ode15s(...
    f, [0 period], x_with_monodromy, integration_opt);
  y_end = trajectory(end,1:nphases)';
  monodromy = trajectory(end,nphases+1:end);
  monodromy = reshape(monodromy, [nphases nphases]);

function dydt_mon = dydt_monodromy_full(t,y, parameters)
  global cds
  y_ode = y(1:cds.nphases);
  
  y_mon = reshape(y(cds.nphases+1:end),cds.nphases,cds.nphases);
  dydt_mon = [
      cds.dydt_ode(t, y_ode, parameters{:}); 
      reshape( ...
        cds.jacobian_ode(t, y_ode, parameters{:}) * y_mon, ...
        [cds.nphases^2 1]) 
  ];

  
  
function [y_end, monodromy] = monodromy_column_by_column(x, period, parameters)
  global cds;
  integration_opt = odeset(...
    'AbsTol',      1e-10,    ...
    'RelTol',      1e-10,    ...
    'BDF',         'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1,      ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  f = @(t, y) cds.dydt_ode(t, y, parameters{:});
  cycle = ode15s(f, linspace(0,period,cds.nDiscretizationPoints), x, integration_opt);
  y_end = deval(cycle,period);
  monodromy = eye(cds.nphases);
  integration_opt = odeset(integration_opt, 'Jacobian', ...
    @(t,y) feval(cds.jacobian_ode, t, deval(cycle,t), parameters{:}));
  f = @(t, y) cds.jacobian_ode(t, deval(cycle,t), parameters{:}) * y;
  parfor i=1:cds.nphases
    fprintf('%d ',i);
    [~, monodromy_map_trajectory] = ode15s(...
      f, [0 period], monodromy(:,i), integration_opt);
    monodromy(:,i) = monodromy_map_trajectory(end,:);
  end 

function init(~,~)
%----------------------------------------------------------
function out=defaultprocessor(varargin)
  global cds
  [y, ~, parameters]           = getComponents(varargin{1}.x);
  cds.previous_phases          = y(1:cds.nphases);
  cds.previous_dydt_0          = cds.dydt_ode(0, ...
                                    cds.previous_phases,parameters{:});
  out                          = varargin{1};
  savePoint(varargin{1});

function options

function CISdata = curve_CIS_first_point(x) %#ok<INUSD>
  CISdata = 1;
  
function CISdata = curve_CIS_step(x, CISdata_in) %#ok<INUSD>
  CISdata = 1;

function adapt(varargin)

function [y,period,parameters] = getComponents(x)
  global cds
  y                            = x(1:cds.nphases*cds.nShootingPoints);
  period                       = x(end-1);
  parameter_value              = x(end);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = parameter_value;
  parameters                   = num2cell(parameters);
  