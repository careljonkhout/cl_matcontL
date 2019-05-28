% continuation of cycles in brusselator
format long
run_init_if_needed
% continuation of cycles cycles in brusselator
odefile = @brusselator_1d; %@brusselator_N_2;
N = 30;
L = 0.5; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
ode_parameters = {N; L; A; B; Dx; Dy};%parameters = {L; A; B; Dx; Dy};
clear global cds
clear global lds
close all
global cds

title_format_string = ...
  'Brusselator N:%d  L:%.2f  A:%.0f  B:%.2f  Dx:%.3f  Dy:%.3f';
title_format_args = {N; L; A; B; Dx; Dy;};

fprintf([title_format_string '\n'], title_format_args{:});



initial_continuation_data = init_multiple_shooting_find_stable_cycle( ...
  'initial_point',             ones(2*N,1), ...
  'time_to_converge_to_cycle', 200, ...
  'subspace_size',             20, ...
  'odefile',                   odefile,  ...
  'ode_parameters',            ode_parameters, ...
  'nMeshIntervals',            2, ...
  'active_parameter_index',    2, ...
  'lower_bound_period',        1, ...
  'upper_bound_period',        20, ...
  'show_plot',                 true ...
);

% 'coordinates_to_plot',       1:30:N ...

disp(initial_continuation_data(end-1))

%% Continue limit cycle from orbit


opt = contset();
opt = contset(opt, 'MaxNumPoints',   10000);
opt = contset(opt, 'InitStepsize',   1e-1);
opt = contset(opt, 'MinStepsize',    1e-10);
opt = contset(opt, 'MaxStepsize',    5e-2);
opt = contset(opt, 'MaxNewtonIters', 3);
opt = contset(opt, 'MaxCorrIters',   6);
opt = contset(opt, 'MaxTestIters',   20);
opt = contset(opt, 'VarTolerance',   1e-6);
opt = contset(opt, 'FunTolerance',   1e-6);
opt = contset(opt, 'NewtonPicardBasisTolerance',   1e-6);
opt = contset(opt, 'contL_SmoothingAngle',   3);
% we don't want to adapt
% since it is not implemented
opt = contset(opt, 'Adapt',          1000*1000*1000);
opt = contset(opt, 'CheckClosed',    1000);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       false);
opt = contset(opt, 'Singularities',  true);
opt = contset(opt, 'CIS_UsingCIS',   false);
opt = contset(opt, 'NewtonPicard',   true);
opt = contset(opt, 'console_output_level',   4);
opt = contset(opt, 'contL_DiagnosticsLevel', 4);
opt = contset(opt, 'PicardTolerance', 1e-8);

[s, datafile] = contL(@multiple_shooting, ...
  initial_continuation_data, ...
  [], opt); 


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
  phases_0                     = xx(1:end-2);
  ode_parameters                   = cds.P0;
  ode_parameters(cds.ActiveParams) = xx(end);
  ode_parameters                   = num2cell(ode_parameters);
  [t,y] = compute_cycle(phases_0,period,ode_parameters);
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

function [value, isterminal, direction] = returnToPlane(t, x)
  global cds;
  % x and should be a column vector
  value = cds.previous_dydt_0'*(x-cds.previous_phases);
  isterminal = t > cds.minimum_period ...
    && max(abs(x-cds.previous_phases)) < cds.poincare_tolerance;
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
