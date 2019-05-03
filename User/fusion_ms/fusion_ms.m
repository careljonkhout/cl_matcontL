
clc
clear global
N = 25;                     
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
a = -1;
b = -0.3;
q_inf = -0.72;

ode_parameters = {a ; b; q_inf};



init_multiple_shooting_find_stable_cycle( ...
  'initial_point',             ones(3*(N-1),1), ...
  'time_to_converge_to_cycle', 150, ... 
  'odefile',                   odefile,  ...
  'ode_parameters',            ode_parameters, ...
  'active_parameter_index',    3, ...
  'nMeshIntervals',            10, ...
  'lower_bound_period',        1, ...
  'upper_bound_period',        20, ...
  'subspace_size',             10, ...
  'show_plot',                 false ...
);

disp(initial_continuation_data(end-1))
opt = contset();
opt = contset(opt, 'MaxNumPoints',            1000);
opt = contset(opt, 'InitStepsize',            0.1);
opt = contset(opt, 'MinStepsize',             1e-6);
opt = contset(opt, 'MaxStepsize',             0.1);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'MaxTestIters',            10);
opt = contset(opt, 'Backward',                false);
opt = contset(opt, 'VarTolerance',            1e-6);
opt = contset(opt, 'FunTolerance',            1e-6);
opt = contset(opt, 'Adapt',                   3);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    5);
opt = contset(opt, 'contL_DiagnosticsLevel',  5);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    100);
opt = contset(opt, 'NewtonPicard',            true);
opt = contset(opt, 'integration_rel_tol',            1e-9);
opt = contset(opt, 'integration_abs_tol',            1e-9);
opt = contset(opt,  ...
                   'every_point_in_separate_mat_file', true, ...
                   'enable_nf_ns',                     false, ...
                   'enable_nf_lpc',                    false, ...
                   'enable_nf_pd',                     false, ...
                   'NewtonPicard',                     true);
 

contL(@multiple_shooting, initial_continuation_data, [], opt);