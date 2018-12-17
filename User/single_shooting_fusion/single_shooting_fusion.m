% continuation of cycles in brusselator

run_init_if_needed
% continuation of cycles cycles in brusselator
N=25;
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));

a = -1; b = -0.3; q_inf = -0.72;
parameters = {a;b;q_inf};
handles = feval(odefile);
title_format_string = 'fusion  N:%d  a:%.1f  b:%.1f  q_{infty}:%.2f';
title_format_args = {N;a;b;q_inf};
clear global cds
global cds
cds.poincare_tolerance = 1e-4;
cds.dydt_ode = handles{2};
cds.jacobian_ode = handles{3};

cds.probfile = odefile;
cds.nap = 1;
cds.ActiveParams = 3;
cds.nphases = 3*(N-1);
cds.ndim = cds.nphases + cds.nap + 1;
cds.P0 = cell2mat(parameters);
cds.options = contset();
cds.options.PartitionMonodromy = true;
cds.nDiscretizationPoints = 1000;
cds.symjac = false;
cds.usernorm = [];
cds.probfile = odefile;
cds.ncoo = cds.nphases;
cds.sout = [];

int_opt = odeset( ...
  'AbsTol',      1e-10,    ...
  'RelTol',      1e-13,    ...
  'Jacobian',     @(t,y) feval(handles{3},t,y,parameters{:}) ...
);


x0 = ones(cds.nphases,1);
f =@(t, y) cds.dydt_ode(t, y, parameters{:});
[t1, x1] = ode15s(f, [0 200], x0, int_opt);

draw_plots = false || true;
if draw_plots
  figure(1)
  plot(t1,x1)
  title(sprintf(title_format_string, title_format_args{:}));
  xlabel('t');
  ylabel('x_1,x_2,y_1,y_2');
end


approximate_period = 12;

cds.x0 = x1(end,:)';
cds.dydt_0 = f(0,x0);
int_opt = odeset(int_opt, 'Events', @returnToPlane);
cds.nDiscretizationPoints = 100;
[t2,x2] = ode15s(f, [0 approximate_period], x1(end,:), int_opt); 
T = t2(end);
int_opt = odeset(int_opt, 'Events', []);
[t3,x3] = ode15s(f, linspace(0,T,cds.nDiscretizationPoints), x2(end,:), int_opt); 

cds.x0_prime = zeros(cds.nDiscretizationPoints,cds.nphases);
for i=1:cds.nDiscretizationPoints
  cds.x0_prime(i,:) = f(0,x3(i,:));
end



% use poincare section to converge to 
%[t3,x3] = ode15s(f, 0:5*approximate_period,x1(end,:)',integration_opt);



if draw_plots
  figure(2)
  plot(t2,x2)
  title(sprintf(title_format_string, title_format_args{:}));
  xlabel('t')
  ylabel('x_1,x_2,y_1,y_2');
end


if draw_plots
  figure(3)
  plot(t3,x3)
  title(sprintf(title_format_string, title_format_args{:}));
  xlabel('t')
  ylabel('x_1,x_2,y_1,y_2');
end


if draw_plots
  figure(4)
  plot(x2(:,1),x2(:,N+1))
  title(sprintf(title_format_string, title_format_args{:}));
  xlabel('x_1')
  ylabel('y_N')
  drawnow
end

%% Continue limit cycle from orbit

opt = contset();
opt = contset(opt, 'InitStepsize',   1e-2);
opt = contset(opt, 'MinStepsize',    1e-4);
opt = contset(opt, 'MaxStepsize',    1e-1);
opt = contset(opt, 'MaxNewtonIters', 3);
opt = contset(opt, 'MaxCorrIters',   10);
opt = contset(opt, 'MaxTestIters',   10);
opt = contset(opt, 'VarTolerance',   1e-6);
opt = contset(opt, 'FunTolerance',   1e-6);
% we don't want to adapt
% since it is not implemented
opt = contset(opt, 'Adapt',          1000*1000*1000);
opt = contset(opt, 'MaxNumPoints',   3);
opt = contset(opt, 'CheckClosed',    50);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       false);
opt = contset(opt, 'Singularities',  false);
opt = contset(opt, 'CIS_UsingCIS',   false);
opt = contset(opt, 'contL_ParallelComputing',   true);


[s, datafile] = contL(@single_shooting,[cds.x0; T; cds.P0(cds.ActiveParams)],[],opt); 


[x,~] = loadPoint(datafile); % DV: load computed cycles
%load('Data\testbruss_Orb_LC.mat')    % DV: load singular points
close all
figure
hold on
for i=1:opt.MaxNumPoints
  xx = x(1:end,i);
  period                       = xx(end-1);
  phases_0                     = xx(1:end-2);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = xx(end);
  parameters                   = num2cell(parameters);
  [t,y] = compute_cycle(phases_0,period,parameters);
  coord1_vals = y(:,1);
  coord2_vals = y(:,3);
  plot(coord1_vals,coord2_vals,'b');
  title(sprintf(title_format_string, title_format_args{:}));
  xlabel('n_1')
  ylabel('U_1')
end


function [t,x] = compute_cycle(x, period, parameters)
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
  [t,x] = ode15s(f, [0 period], x, integration_opt);
end

function [value, isterminal, direction] = returnToPlane(t, x)
  global cds;
  % x and should be a column vectors
  value = cds.dydt_0'*(x-cds.x0);
  isterminal = t > 1 && sum((x-cds.x0).^2) < cds.poincare_tolerance;
  direction = 1;
end
