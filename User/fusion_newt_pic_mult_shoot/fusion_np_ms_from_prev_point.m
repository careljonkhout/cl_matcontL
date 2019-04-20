% continuation of cycles in fusion system


run_init_if_needed

N = 25;
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
handles = feval(odefile);
a = -1; b = -0.3; q_inf = -0.72;
parameters = {a;b;q_inf};
title_string_args = {N,a,b,q_inf};

clear global cds
clear global lds

global cds
cds.poincare_tolerance = 1e-4;
cds.minimum_period = 1;
cds.nphases = 3*(N-1);
cds.nap = 1;
cds.ActiveParams = 3;
cds.ndim = cds.nphases + cds.nap + 1;
cds.P0 = cell2mat(parameters);
cds.options = contset();
cds.options.PartitionMonodromy = true;
cds.preferred_basis_size = 5;
cds.symjac = true;
cds.usernorm = [];
cds.probfile = odefile;
cds.dydt_ode = handles{2};
cds.jacobian_ode = handles{3};
cds.ncoo = cds.nphases;


a = -1;
b = -0.3;
q_inf = -0.72;
parameters = {a;b;q_inf};

load('point 16.mat');
parameters{3} = point.x(end);



int_opt = odeset( ...
  'Jacobian',     @(t,y) feval(handles{3},t,y,parameters{:}) ...
);


x0 = ones(cds.nphases,1);
dydt = handles{2};
f =@(t, y) dydt(t, y, parameters{:});

[t1, x1] = ode15s(f, [0 400], point.x(1:end-2), int_opt);


draw_plots = false || true;
if draw_plots
  figure(1)
  plot(t1,x1)
  title('fusion');
  xlabel('t');
  ylabel('values at grid points');
end





opt = contset();
opt = contset(opt, 'InitStepsize',   1e-1);
opt = contset(opt, 'MinStepsize',    1e-10);
opt = contset(opt, 'MaxStepsize',    1e-1);
opt = contset(opt, 'MaxNewtonIters', 8);
opt = contset(opt, 'MaxCorrIters',   10);
opt = contset(opt, 'MaxTestIters',   10);
opt = contset(opt, 'VarTolerance',   1e-6);
opt = contset(opt, 'FunTolerance',   1e-6);
% we don't want to adapt
% since it is not implemented
opt = contset(opt, 'Adapt',          1000*1000*1000);
opt = contset(opt, 'MaxNumPoints',   100*1000);
opt = contset(opt, 'contL_SmoothingAngle', 10);
opt = contset(opt, 'CheckClosed',    50000);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       false);
opt = contset(opt, 'Singularities',  true);
opt = contset(opt, 'CIS_UsingCIS',   false);
opt = contset(opt, 'NewtonPicard',   true);
opt = contset(opt, 'every_point_in_separate_mat_file', true);
opt = contset(opt, 'console_output_level',   5);
opt = contset(opt, 'contL_DiagnosticsLevel', 5);
opt = contset(opt, 'Multipliers', true);

tol = 1e-9;

opt = contset(opt, 'integration_abs_tol', tol);
opt = contset(opt, 'integration_rel_tol', tol);
opt = contset(opt, 'multipliers_abs_tol', 1e-10);
opt = contset(opt, 'multipliers_rel_tol', 1e-10);

orbit                   = x1;
nMeshPoints         = 20;
active_parameter_index  = 3;
time_integration_method = @ode15s;
lower_bound_period      = 1;
upper_bound_period      = 30;

time_integration_options = odeset( ...
    'AbsTol',      tol,    ...
    'RelTol',      tol);

initial_continuation_data = init_multiple_shooting_from_orbit(...
    orbit, ...
    nMeshPoints, ...
    odefile, ...
    parameters, ...
    active_parameter_index, ...
    time_integration_method, ...
    lower_bound_period, ...
    upper_bound_period, ...
    time_integration_options);
[s, datafile] = contL(@multiple_shooting, ...
  initial_continuation_data,[], opt); 



% use poincare section to converge to 
%[t3,x3] = ode15s(f, 0:5*approximate_period,x1(end,:)',integration_opt);


if draw_plots
  figure(2)
  hold on;
  plot(t2,x2)
   title(sprintf( ...
    'fusion N:%d a:%.2f b:%.2f q_{inf}:%.2f', ...
    title_string_args{:}));
  xlabel('t')
  ylabel('x_1,x_2,y_1,y_2');
end

if draw_plots
  figure(3)
  hold on;
  plot(t3,x3)
   title(sprintf( ...
    'fusion N:%d a:%.2f b:%.2f q_{inf}:%.2f', ...
    title_string_args{:}));
  xlabel('t')
  ylabel('n_1, ..., n_{N-1},U_1, ..., U_{N-1},z_1,...,z_{N-1}');
end


if draw_plots
  figure(3)
  plot(x2(:,1),x2(:,2))
  title(sprintf( ...
    'fusion N:%d a:%.2f b:%.2f q_{inf}:%.2f', ...
    title_string_args{:}));
  xlabel('n_1')
  ylabel('U_1')

  drawnow
end

%% Continue limit cycle from orbit




[s, datafile] = contL(@multiple_shooting, ...
  [cds.previous_phases; period; cds.P0(cds.ActiveParams)],[],opt); 


figure
hold on;
coordinate1 = 1;
coordinate2 = 2;
 title(sprintf( ...
    'fusion N:%d a:%.2f n:%.2f q_{inf}:%.2f', ...
     N,a,b,q_inf));
 xlabel('x_1')
 ylabel('y_1')
 [x, v, h, mult] = loadPoint(datafile); % DV: load computed cycles
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


function [value, isterminal, direction] = returnToPlane(t, x)
  global cds;
  % x and should be a column vector
  value = cds.previous_dydt_0'*(x-cds.previous_phases);
  isterminal = t > cds.minimum_period ...
    && max(abs(x-cds.previous_phases)) < cds.poincare_tolerance;
  direction = 1;
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

