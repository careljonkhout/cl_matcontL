% test of Neimark Sacker bifurcation detection using orhtogonal collocation
clc
clear global
N = 300;                     
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d',N));
Gamma_inf = -0.8;
q_inf     = -0.72;
D0        = 1.9;
D1        = -1.1;
D2        = 0;
a         = -1;
b         = -0.3;
zeta1     = 1.1;
mu1       = 0.05;
epsilon   = 0.05;
ZS        = 0;
gamma1    = 5/3;
lambda_n  = 1.25;
lambda_T  = 1.5;
c_n       = 1.1;
c_T       = 0.9;

%ode_parameters ={ 
%  N, Gamma_inf, q_inf, D0, D1, D2, a, b, zeta1, mu1, epsilon, ZS, gamma1, ...
%  lambda_n, lambda_T, c_n, c_T};
ode_parameters = { a,b,q_inf};




initial_continuation_data = init_single_shooting_stable_cycle( ...
  'initial_point',             ones(3*(N-1),1), ...
  'time_to_converge_to_cycle', 150, ... 
  'odefile',                   odefile,  ...
  'ode_parameters',            ode_parameters, ...
  'active_parameter_index',    3, ...
  'lower_bound_period',        1, ...
  'upper_bound_period',        20, ...
  'subspace_size',             10, ...
  'show_plot',                 false ...
);

disp(initial_continuation_data(end-1))
opt = contset();
opt = contset(opt, 'MaxNumPoints',            100);
opt = contset(opt, 'InitStepsize',            0.1);
opt = contset(opt, 'MinStepsize',             1e-6);
opt = contset(opt, 'MaxStepsize',             0.1);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'MaxTestIters',            10);
opt = contset(opt, 'Backward',                false);
opt = contset(opt, 'VarTolerance',            1e-3);
opt = contset(opt, 'FunTolerance',            1e-3);
opt = contset(opt, 'Adapt',                   3);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    5);
opt = contset(opt, 'contL_DiagnosticsLevel',  5);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    1);
opt = contset(opt, 'NewtonPicard',            true);
opt = contset(opt, 'integration_rel_tol',            1e-6);
opt = contset(opt, 'integration_abs_tol',            1e-6);
opt = contset(opt, 'enable_bpc',              false);
 

contL(@single_shooting, initial_continuation_data, [], opt);