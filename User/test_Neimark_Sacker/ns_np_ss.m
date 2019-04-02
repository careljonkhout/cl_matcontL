% test of Neimark Sacker bifurcation detection
clc
clear global
load('H_LC.mat')
step = 20;
ode_parameters                           = num2cell(lds.P0);
ode_parameters{lds.ActiveParams} = x(end,step);

initial_continuation_data = init_single_shooting( ...
  'point_on_limitcycle',     x(1:4, step), ...
  'odefile',                 @STLARTEST,  ...
  'ode_parameters',          ode_parameters, ...
  'active_parameter_index',  lds.ActiveParams, ...
  'lower_bound_period',      1, ...
  'upper_bound_period',      x(end-1, step)+1, ...
  'subspace_size',           4);

opt = contset();
opt = contset(opt, 'MaxNumPoints',   30);
opt = contset(opt, 'InitStepsize',   1e-1);
opt = contset(opt, 'MinStepsize',    1e-10);
opt = contset(opt, 'MaxStepsize',    1);
opt = contset(opt, 'MaxNewtonIters', 8);
opt = contset(opt, 'MaxCorrIters',   10);
opt = contset(opt, 'MaxTestIters',   10);
opt = contset(opt, 'VarTolerance',   1e-5);
opt = contset(opt, 'FunTolerance',   1e-9);
opt = contset(opt, 'Adapt',          1000*1000*1000);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Singularities',          true);
opt = contset(opt, 'console_output_level',   5);
opt = contset(opt, 'contL_DiagnosticsLevel', 5);
opt = contset(opt, 'MoorePenrose'           , false);
opt = contset(opt, 'integration_abs_tol'        , 1e-13);
opt = contset(opt, 'integration_rel_tol'        , 1e-13);




global contopts
contopts = opt;

ss_handles = single_shooting();
feval(ss_handles{1},initial_continuation_data)
 

contL(@single_shooting, initial_continuation_data, [], opt);