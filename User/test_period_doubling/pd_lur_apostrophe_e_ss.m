% test of period doubling bifurcation (PD) detection using single_shooting
% also tests the locating routines for period doubling
clc
clear global
load('H_LC.mat')

step = 10;
ode_parameters                   = num2cell(lds.P0);
ode_parameters{lds.ActiveParams} = x(end,step);
disp(x(end-1,step))

initial_continuation_data = init_single_shooting( ...
  'point_on_limitcycle',     x(1:3,step), ...
  'odefile',                 @Lur_apostrophe_e,  ...
  'ode_parameters',          ode_parameters, ...
  'active_parameter_index',  lds.ActiveParams, ...
  'lower_bound_period',      x(end-1,step)-1, ...
  'upper_bound_period',      x(end-1,step)+1, ...
  'subspace_size',           3, ...
  'show_plot',               false ...
);

disp(initial_continuation_data(end-1))
opt = contset();
opt = contset(opt, 'MaxNumPoints',            25);
opt = contset(opt, 'InitStepsize',            1e-1);
opt = contset(opt, 'MinStepsize',             1e-10);
opt = contset(opt, 'MaxStepsize',             1e-1);
opt = contset(opt, 'Backward',                true);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'VarTolerance',            1e-5);
opt = contset(opt, 'FunTolerance',            1e-9);
opt = contset(opt, 'MaxTestIters',            20);
opt = contset(opt, 'contL_Testf_VarTolerance', 1e-6);
opt = contset(opt, 'contL_Testf_FunTolerance', 1e-6);
opt = contset(opt, 'Adapt',                   1000*1000*1000);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    3);
opt = contset(opt, 'contL_DiagnosticsLevel',  0);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    pi);

global contopts
contopts = opt;

% Period Doubling (period = 8.386067e+00, parameter = 6.273247e-01)
contL(@single_shooting, initial_continuation_data, [], opt);