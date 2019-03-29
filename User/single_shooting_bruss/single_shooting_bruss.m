% continuation of cycles in brusselator
clear global cds
clear global lds
run_init_if_needed
% continuation of cycles cycles in brusselator
odefile = @brusselator_N_2;
N=2;
L = 1.1; A = 1; B = 2.2; Dx = 0.008; Dy = 0.004;
parameters = {L; A; B; Dx; Dy};
handles = feval(odefile);
title_format_string = ...
  'Brusselator N:%d  L:%.0f  A:%.0f  B:%.1f  Dx:%.3f  Dy:%.3f';
title_format_args = {N; L; A; B; Dx; Dy;};
clear global cds
clear global lds
global cds
cds.poincare_tolerance = 1e-4;
cds.dydt_ode = handles{2};
cds.jacobian_ode = handles{3};
cds.jacobian_params = handles{4};
cds.p = 4;
cds.integrator = @ode15s;
cds.probfile = odefile;
cds.nap = 1;
cds.ActiveParams = 4;
cds.nphases = 2*N;
cds.ndim = cds.nphases + cds.nap + 1;
cds.P0 = cell2mat(parameters);
cds.options = contset();
cds.options.PartitionMonodromy = false;
cds.nDiscretizationPoints = 1000;
cds.symjac = false;
cds.usernorm = [];
cds.probfile = odefile;
cds.ncoo = cds.nphases;
cds.mv_count = 0;

    
int_opt = odeset( ...
  'AbsTol',      1e-10,    ...
  'RelTol',      1e-13,    ...
  'Jacobian',     @(t,y) feval(handles{3},t,y,parameters{:}) ...
);


x0 = ones(cds.nphases,1);
dydt = handles{2};
f =@(t, y) dydt(t, y, parameters{:});
[t1, x1] = ode15s(f, [0 150], x0, int_opt);

draw_plots = false || false;
if draw_plots
  figure(1)
  plot(t1,x1)
  title(sprintf(title_format_string, title_format_args{:}));
  xlabel('t');
  ylabel('x_1,x_2,y_1,y_2');
end

approximate_period = 12;

cds.previous_phases = x1(end,:)';
cds.previous_dydt_0 = f(0,x0);
int_opt = odeset(int_opt, 'Events', @returnToPlane);

[t2,x2] = ode15s(f, linspace(0,approximate_period,1000), x1(end,:), int_opt); 
period = t2(end);
int_opt = odeset(int_opt, 'Events', []);

if draw_plots
  figure(2)
  plot(t2,x2)
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
opt = contset(opt, 'MaxNumPoints',   17);
opt = contset(opt, 'InitStepsize',   1e-2);
opt = contset(opt, 'MinStepsize',    1e-6);
opt = contset(opt, 'MaxStepsize',    1e-1);
opt = contset(opt, 'MaxNewtonIters', 3);
opt = contset(opt, 'MaxCorrIters',   10);
opt = contset(opt, 'MaxTestIters',   10);
opt = contset(opt, 'VarTolerance',   1e-6);
opt = contset(opt, 'FunTolerance',   1e-6);
% we don't want to adapt
% since it is not implemented
opt = contset(opt, 'Adapt',          1000*1000*1000);
opt = contset(opt, 'CheckClosed',    50);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       false);
opt = contset(opt, 'Singularities',  false);
opt = contset(opt, 'CIS_UsingCIS',   false);
opt = contset(opt, 'contL_SmoothingAngle',   3);
opt = contset(opt, 'console_output_level',   5);
opt = contset(opt, 'contL_DiagnosticsLevel', 5);

initial_continuation_data = [x0; period; cds.P0(cds.ActiveParams)];
initial_continuation_tangent_vector = [];
[s, datafile] = contL(@single_shooting, ...
  initial_continuation_data, ...
  initial_continuation_tangent_vector, opt); 


figure
hold on;
coordinate1 = 1;
coordinate2 = 2;
 title(sprintf( ...
    'brusselator'))
 ylabel('y_1')
x = loadPoint(datafile); % DV: load computed cycles
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

return
figure
plot(x(end,:),x(end-1,:))
title(sprintf(title_format_string, title_format_args{:}));
xlabel('TODO: set xlabel')
ylabel('period')

figure
hold on;
nMults = size(mult,1);
for i=nMults-10:nMults
  plot(x(end,:),mult(i, :))
end

function [value, isterminal, direction] = returnToPlane(t, x)
  global cds;
  % x and should be a column vectors
  value = cds.previous_dydt_0'*(x-cds.previous_phases);
  isterminal = t > 1 ...
    && sum((x-cds.previous_phases).^2) < cds.poincare_tolerance;
  direction = 1;
end

function [time_values, trajectory] ...
  = compute_cycle(initial_condition, period, parameters)
  global cds
  f =@(t, y) cds.dydt_ode(t, y, parameters{:});
  integration_opt = odeset(...
    'AbsTol',      1e-10,    ...
    'RelTol',      1e-10,    ...
    'BDF',         'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1,      ...
    'Events',      [],...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  time_values = linspace(0,period, 1000);
  [time_values, trajectory] ...
    = ode15s(f, time_values, initial_condition, integration_opt);
end