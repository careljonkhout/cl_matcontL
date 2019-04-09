% test of Neimark Sacker bifurcation detection using orhtogonal collocation
clc
clear global
N = 50;                     
odefile = @nonadiabatic_tubular_reactor;
D     = 0.165;
P_em  = 5;
P_eh  = 5;
BETA  = 2.35;
phi_0 = 1;
GAMMA = 25;
B     = 0.5;

ode_parameters = {D; P_em; P_eh; BETA; phi_0; GAMMA; B};



initial_continuation_data = initOrbLC_L_find_stable_cycle( ...
  'initial_point',             ones(2*N,1), ...
  'time_to_converge_to_cycle', 150, ... 
  'odefile',                   odefile,  ...
  'ode_parameters',            ode_parameters, ...
  'active_parameter_index',    1, ...
  'lower_bound_period',        1, ...
  'upper_bound_period',        20, ...
  'nMeshIntervals',            40, ...
  'show_plot',                 false ...
);

disp(initial_continuation_data(end-1))
opt = contset();
opt = contset(opt, 'MaxNumPoints',            100);
opt = contset(opt, 'InitStepsize',            0.5);
opt = contset(opt, 'MinStepsize',             1e-6);
opt = contset(opt, 'MaxStepsize',             0.5);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'MaxTestIters',            10);
opt = contset(opt, 'Backward',                false);
opt = contset(opt, 'VarTolerance',            1e-6);
opt = contset(opt, 'FunTolerance',            1e-6);
opt = contset(opt, 'Adapt',                   3);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    0);
opt = contset(opt, 'contL_DiagnosticsLevel',  0);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    1);
opt = contset(opt, 'enable_bpc',              false);
opt = contset(opt, 'enable_nf_lpc',           false, ...
                   'every_point_in_separate_mat_file', true);

 

contL(@limitcycleL, initial_continuation_data, [], opt);