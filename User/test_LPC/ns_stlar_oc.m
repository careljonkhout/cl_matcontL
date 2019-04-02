% test of Neimark Sacker bifurcation detection using orhtogonal collocation
clc
clear global
load('H_LC.mat')

step = size(x,2);
ode_parameters                           = num2cell(lds.P0);
ode_parameters{lds.ActiveParams} = x(end,step);
disp(x(end-1,step))


initial_continuation_data = initOrbLC_L_key_value( ...
  'point_on_limitcycle',     x(1:4,step), ...
  'odefile',                 @STLARTEST,  ...
  'ode_parameters',          ode_parameters, ...
  'active_parameter_index',  lds.ActiveParams, ...
  'lower_bound_period',      1, ...
  'upper_bound_period',      x(end-1,step)+1, ...
  'nMeshIntervals',          100, ...
  'show_plot',               false ...
);

disp(initial_continuation_data(end-1))
opt = contset();
opt = contset(opt, 'MaxNumPoints',            100);
opt = contset(opt, 'InitStepsize',            1);
opt = contset(opt, 'MinStepsize',             1e-6);
opt = contset(opt, 'MaxStepsize',             1);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'MaxTestIters',            10);
opt = contset(opt, 'Backward',                false);
opt = contset(opt, 'VarTolerance',            1e-3);
opt = contset(opt, 'FunTolerance',            1e-3);
opt = contset(opt, 'Adapt',                   3);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    0);
opt = contset(opt, 'contL_DiagnosticsLevel',  0);
opt = contset(opt, 'MoorePenrose',            false);
%opt = contset(opt, 'contL_SmoothingAngle',    0.1);

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