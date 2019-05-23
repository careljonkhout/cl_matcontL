% test of limit point bifurcation detection using orhtogonal collocation
clc
clear global
N = 50;                     
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
a = -1;
b = -0.3;  
q_inf = -0.72;

ode_parameters = {a ; b; q_inf};


initial_continuation_data = init_collocation_find_stable_cycle( ...
  'initial_point',             ones(3*(N-1),1), ...
  'time_to_converge_to_cycle', 150, ... 
  'odefile',                   odefile,  ...
  'ode_parameters',            ode_parameters, ...
  'active_parameter_index',    3, ...
  'lower_bound_period',        1, ...
  'upper_bound_period',        20, ...
  'nMeshIntervals',            40, ...
  'show_plot',                 false ...
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
opt = contset(opt, 'VarTolerance',            1e-3);
opt = contset(opt, 'FunTolerance',            1e-3);
opt = contset(opt, 'Adapt',                   3);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    3);
opt = contset(opt, 'contL_DiagnosticsLevel',  3);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    1);
opt = contset(opt, ...
                   'enable_nf_pd',            false, ...
                   'enable_nf_lpc',           true, ...
                   'enable_nf_ns',            false);

 

contL(@limitcycleL, initial_continuation_data, [], opt, 'callback', @plot_T_versus_param);