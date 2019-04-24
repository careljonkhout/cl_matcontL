

clc
clear global
N = 75;                     
%odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
odefile = str2func(sprintf('fusion_mex_N_%d', N));
a = -1;
b = -0.3;  
q_inf = -0.72;

ode_parameters = {a ; b; q_inf};

initial_continuation_data = init_multiple_shooting_find_stable_cycle( ...
  'initial_point',             ones(3*(N-1),1), ...
  'time_to_converge_to_cycle', 150, ... 
  'subspace_size',             8, ...
  'odefile',                   odefile,  ...
  'ode_parameters',            ode_parameters, ...
  'active_parameter_index',    3, ...
  'lower_bound_period',        1, ...
  'upper_bound_period',        6, ...
  'nMeshIntervals',            20, ...
  'show_plot',                 false ...
);

disp(initial_continuation_data(end-1))
opt = contset();
opt = contset(opt, 'MaxNumPoints',            1000);
opt = contset(opt, 'InitStepsize',            0.25);
opt = contset(opt, 'MinStepsize',             1e-6);
opt = contset(opt, 'MaxStepsize',             0.25);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'MaxTestIters',            10);
opt = contset(opt, 'Backward',                false);
opt = contset(opt, 'VarTolerance',            1e-3);
opt = contset(opt, 'FunTolerance',            1e-3);
opt = contset(opt, 'Adapt',                   3);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'NewtonPicard',            true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    5);
opt = contset(opt, 'contL_DiagnosticsLevel',  5);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    1);
opt = contset(opt, 'enable_bpc',              false, ...
                   'enable_nf_pd',            false, ...
                   'enable_nf_lpc',           true, ...
                   'enable_nf_ns',            false, ...
                   'integration_rel_tol',     1e-8, ...
                   'integration_abs_tol',     1e-8);

 
tic
contL(@multiple_shooting, ...
      initial_continuation_data, [], opt, @plot_T_versus_param);
toc