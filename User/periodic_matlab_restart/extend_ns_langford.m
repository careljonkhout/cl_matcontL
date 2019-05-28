% test of Neimark Sacker bifurcation detection using single_shooting
clc
clear global
dirname  = fullfile(get_path(),'Data','langford_ss');
filename = get_latest_file(dirname);
load(fullfile(dirname,filename));


initial_continuation_data = init_single_shooting( ...
  'point_on_limitcycle',     point.x(1:end-2), ...
  'odefile',                 @Langford,  ...
  'ode_parameters',          ode_parameters, ...
  'active_parameter_index',  globals.lds.ActiveParams, ...
  'lower_bound_period',      1, ...
  'upper_bound_period',      x(end-1,step)+1, ...
  'subspace_size',           3, ...
  'show_plot',               false ...
);

disp(initial_continuation_data(end-1))
opt = contset();
opt = contset(opt, 'MaxNumPoints',            1);
opt = contset(opt, 'InitStepsize',            4e-2);
opt = contset(opt, 'MinStepsize',             1e-10);
opt = contset(opt, 'MaxStepsize',             4e-2);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'VarTolerance',            1e-5);
opt = contset(opt, 'FunTolerance',            1e-9);
opt = contset(opt, 'MaxTestIters',            40);
opt = contset(opt, 'contL_Testf_VarTolerance', 1e-6);
opt = contset(opt, 'contL_Testf_FunTolerance', 1e-6);
opt = contset(opt, 'Adapt',                   1000*1000*1000);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    5);
opt = contset(opt, 'contL_DiagnosticsLevel',  0);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    pi, ...
                   'pause',                   true, ...
                   'nsteps_before_pause',     20, ...
                   'Filename',                'langford_ss');

global contopts
contopts = opt;


contL(@single_shooting, initial_continuation_data, [], opt);