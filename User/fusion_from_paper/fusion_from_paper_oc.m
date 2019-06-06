clc
clear global
N = 50;                     
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
a = -1;
b = -0.3;  
q_inf = -0.7;

ode_parameters = {a ; b; q_inf};
 

initial_continuation_data = init_collocation_find_stable_cycle( ...
  'initial_point',             ones(3*(N-1),1), ...
  'time_to_converge_to_cycle', 500, ... 
  'odefile',                   odefile,  ...
  'ode_parameters',            ode_parameters, ...
  'active_parameter_index',    3, ...
  'lower_bound_period',        1, ...
  'upper_bound_period',        70, ...
  'nMeshIntervals',            40, ...
  'nCollocationPoints',        4, ...
  'show_plot',                 true, ...
  'poincare_tolerance',        8e-3 ...
);



fprintf('period: %.4f\n', initial_continuation_data(end-1));

opt = contset( ...
    'MaxNumPoints',                     100000, ...
    'InitStepsize',                     2e-2, ... 
    'MinStepsize',                      1e-8, ...
    'MaxStepsize',                      5, ...
    'MaxNewtonIters',                   8, ...
    'MaxCorrIters',                     10, ...
    'MaxTestIters',                     10, ...
    'Backward',                         false, ...
    'VarTolerance',                     1e-6, ...
    'FunTolerance',                     1e-6, ...
    'Adapt',                            3, ...
    'Multipliers',                      true, ...
    'Singularities',                    true, ...
    'console_output_level',             5, ...
    'contL_DiagnosticsLevel',           5, ...
    'MoorePenrose',                     false, ...
    'newtcorrL_use_max_norm',           true, ...
    'contL_SmoothingAngle',             6 * pi / 180, ...
    'multiplier_print_threshold',       0.95, ...
    'enable_nf_pd',                     false, ...
    'enable_nf_lpc',                    false, ...
    'enable_nf_ns',                     false, ...
    'MoorePenrose',                     false);

 
figure
global cds
cds.usernorm = @(x) max(abs(x));

contL(@limitcycleL, initial_continuation_data, [], opt, ...
            'callback', @plot_T_versus_param);