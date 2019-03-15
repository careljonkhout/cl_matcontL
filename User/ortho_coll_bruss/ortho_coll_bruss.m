% continuation of cycles in brusselator
format long
run_init_if_needed
% continuation of cycles cycles in brusselator
odefile = @brusselator_1d_no_jac; %@brusselator_N_2;
N=200;
L = 1.1; A = 1; B = 2.2; Dx = 0.008; Dy = 0.004;
parameters = {N; L; A; B; Dx; Dy};%parameters = {L; A; B; Dx; Dy};
clear global cds
clear global lds
global cds
handles = feval(odefile);
cds.dydt_ode = handles{2};
active_parameter_index = 2;
handles = feval(odefile);
title_format_string = ...
  'Brusselator N:%d  L:%.0f  A:%.0f  B:%.1f  Dx:%.3f  Dy:%.3f';
title_format_args = {N; L; A; B; Dx; Dy;};


nphases = 2*N;


    
int_opt = odeset( ...
  'AbsTol',      1e-10,    ...
  'RelTol',      1e-13 ...
);
%'Jacobian',     @(t,y) feval(handles{3},t,y,parameters{:}) ...

x0 = ones(nphases,1);
dydt = handles{2};
f =@(t, y) dydt(t, y, parameters{:});
[t1, x1] = ode15s(f, [0 150], x0, int_opt);

draw_plots = true;
if draw_plots
  figure(1)
  plot(t1,x1)
  title(sprintf(title_format_string, title_format_args{:}));
  xlabel('t');
  ylabel('x_1,x_2,y_1,y_2');
end


% approximate_period should between the period and twice the period for the
% orthogonal collocation initializer

approximate_period = 7; 


[t2,x2] = ode15s(f, linspace(0,approximate_period,500) , x1(end,:), int_opt); 

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
opt = contset(opt, 'MaxNumPoints',   4);
opt = contset(opt, 'InitStepsize',   10);
opt = contset(opt, 'MinStepsize',    1e-6);
opt = contset(opt, 'MaxStepsize',    10);
opt = contset(opt, 'MaxNewtonIters', 3);
opt = contset(opt, 'MaxCorrIters',   4);
opt = contset(opt, 'MaxTestIters',   30);
opt = contset(opt, 'VarTolerance',   1e-6);
opt = contset(opt, 'FunTolerance',   1e-6);
opt = contset(opt, 'contL_SmoothingAngle',   3);
% we don't want to adapt
% since it is not implemented
opt = contset(opt, 'Adapt',          1000*1000*1000);
opt = contset(opt, 'CheckClosed',    1000);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       true);
opt = contset(opt, 'Singularities',  false);
opt = contset(opt, 'CIS_UsingCIS',   false);
opt = contset(opt, 'console_output_level',   5);

ntst = 20; % number of mesh intervals
ncol = 4;  % number of collocation points
init_orb_lc_tolerance = 1e-2;

[x0_for_contL,~] = initOrbLC(odefile, t2, x2, cell2mat(parameters), ...
  active_parameter_index, ntst, ncol, init_orb_lc_tolerance);

[sout, datafile] = contL(@limitcycleL,x0_for_contL,[],opt); %#ok<*ASGLU>



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
for i=1:size(x,2)
  xx = x(1:end,i);
  period                       = xx(end-1);
  phases_0                     = xx(1:2*N);
  parameters                   = cds.P0;
  parameters(cds.ActiveParams) = xx(end);
  parameters                   = num2cell(parameters);
  [t,y] = compute_cycle(phases_0,period,parameters);
  coord1_vals = y(:,1);
  coord2_vals = y(:,end);
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
    'Events',      [] ...
  );
  time_values = linspace(0,period, 1000);
  [time_values, trajectory] ...
    = ode15s(f, time_values, initial_condition, integration_opt);
end