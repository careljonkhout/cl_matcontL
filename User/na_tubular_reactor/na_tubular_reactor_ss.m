clc
clear global
close all;
N = 20;                     
odefile = @nonadiabatic_tubular_reactor;
D     = 0.156;
P_em  = 5;
P_eh  = 5;
BETA  = 2.35;
phi_0 = 1;
GAMMA = 25;
B     = 0.5;

ode_parameters = {D; P_em; P_eh; BETA; phi_0; GAMMA; B};



initial_continuation_data = init_single_shooting_stable_cycle( ...
  'initial_point',             ones(2*N,1), ...
  'time_to_converge_to_cycle', 150, ... 
  'subspace_size',             10, ...
  'odefile',                   odefile,  ...
  'ode_parameters',            ode_parameters, ...
  'active_parameter_index',    1, ...
  'lower_bound_period',        1, ...
  'upper_bound_period',        20, ...
  'show_plot',                 false ...
);

disp(initial_continuation_data(end-1))
opt = contset();
opt = contset(opt, 'MaxNumPoints',            4000);
opt = contset(opt, 'InitStepsize',            0.3);
opt = contset(opt, 'MinStepsize',             1e-6);
opt = contset(opt, 'MaxStepsize',             0.3);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'MaxTestIters',            10);
opt = contset(opt, 'Backward',                true);
opt = contset(opt, 'VarTolerance',            1e-6);
opt = contset(opt, 'FunTolerance',            1e-6);
opt = contset(opt, 'Adapt',                   3);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    3);
opt = contset(opt, 'contL_DiagnosticsLevel',  3);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    1);
opt = contset(opt, 'NewtonPicard',            true);

opt = contset(opt, 'enable_nf_lpc',           true, ...
                   'enable_nf_pd',            false, ...
                   'enable_nf_ns',            false);
opt = contset(opt, 'every_point_in_separate_mat_file', true);

 

contL(@single_shooting, initial_continuation_data, [], opt, @plot_T_versus_param);