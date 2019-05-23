% test of limit point bifurcation detection using orhtogonal collocation
clc
clear global
N = 50;                     
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
a = -1;
b = -0.3;  
q_inf = -0.7; 

dirname = [mfilename '_N_' num2str(N)];

ode_parameters = {a ; b; q_inf};


initial_continuation_data = init_single_shooting_find_stable_cycle( ...
  'initial_point',             ones(3*(N-1),1), ...
  'time_to_converge_to_cycle', 1000, ... 
  'odefile',                   odefile,  ...
  'ode_parameters',            ode_parameters, ...
  'subspace_size',             15, ...
  'active_parameter_index',    3, ...
  'lower_bound_period',        1, ...
  'upper_bound_period',        70, ...
  'show_plot',                 true, ...
  'poincare_tolerance',        8e-4, ...
  'time_integration_options',   odeset('AbsTol', 1e-12, 'RelTol', 1e-12) ...
);

disp(initial_continuation_data(end-1))
opt = contset();
opt = contset(opt, 'MaxNumPoints',            5000);
opt = contset(opt, 'InitStepsize',            0.25); 
opt = contset(opt, 'MinStepsize',             1e-6);
opt = contset(opt, 'MaxStepsize',             0.25);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'MaxTestIters',            10);
opt = contset(opt, 'Backward',                false);
opt = contset(opt, 'VarTolerance',            1e-6);
opt = contset(opt, 'FunTolerance',            1e-6);
opt = contset(opt, 'Adapt',                   3);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    3);
opt = contset(opt, 'contL_DiagnosticsLevel',  3);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    1);
opt = contset(opt,  ...
                   'multiplier_print_threshold', 0.8, ...
                   'enable_nf_pd',            false, ...
                   'enable_nf_lpc',           true, ...
                   'enable_nf_ns',            false, ...
                   'every_point_in_separate_mat_file', true, ...
                   'NewtonPicard',                     true, ...
                   'integration_abs_tol',        1e-12, ...
                   'integration_rel_tol',        1e-12);

 

contL(@single_shooting, initial_continuation_data, [], opt);