% test of Neimark Sacker bifurcation detection using orhtogonal collocation
clc
clear global
load('LC_LC_langford.mat')

step = 1;
ode_parameters                   = num2cell(globals.lds.P0);
ode_parameters{globals.lds.ActiveParams} = x(end,step);
disp(x(end-1,step))


initial_continuation_data = initOrbLC_L_key_value( ...
  'point_on_limitcycle',     x(1:3,step), ...
  'odefile',                 @Langford,  ...
  'ode_parameters',          ode_parameters, ...
  'active_parameter_index',  globals.lds.ActiveParams, ...
  'lower_bound_period',      1, ...
  'upper_bound_period',      x(end-1,step)+1, ...
  'nMeshIntervals',          100, ...
  'show_plot',               false ...
);

disp(initial_continuation_data(end-1))
opt = contset();
opt = contset(opt, 'MaxNumPoints',            30);
opt = contset(opt, 'InitStepsize',            1e-1);
opt = contset(opt, 'MinStepsize',             1e-10);
opt = contset(opt, 'MaxStepsize',             1);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'VarTolerance',            1e-5);
opt = contset(opt, 'FunTolerance',            1e-9);
opt = contset(opt, 'MaxTestIters',            40);
opt = contset(opt, 'contL_Testf_VarTolerance', 1e-9);
opt = contset(opt, 'contL_Testf_FunTolerance', 1e-9);
opt = contset(opt, 'Adapt',                   1000*1000*1000);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    3);
opt = contset(opt, 'contL_DiagnosticsLevel',  3);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    pi);

options.orbit_abs_tol             = 1e-12;
options.orbit_rel_tol             = 1e-12;
options.MV_abs_tol                = 1e-11;
options.MV_rel_tol                = 1e-11;
options.shoot_abs_tol             = 1e-11;
options.shoot_rel_tol             = 1e-11;
options.monodromy_map_abs_tol     = 1e-11;
options.monodromy_map_rel_tol     = 1e-11;
options.continue_subspace_abs_tol = 1e-12;
options.continue_subspace_rel_tol = 1e-12;
options.jacobian_abs_tol          = 1e-12;
options.jacobian_rel_tol          = 1e-12;
options.int_abs_tol               = 1e-9;
options.int_rel_tol               = 1e-9;

global contopts
contopts = opt;


contL(@limitcycleL, initial_continuation_data, [], opt);