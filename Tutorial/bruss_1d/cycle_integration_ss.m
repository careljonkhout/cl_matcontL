% continuation of cycles in brusselator using multiple shooting
format long
odefile = @brusselator_1d;
N = 30;
L = 0.5; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
ode_parameters = {N; L; A; B; Dx; Dy};
clear global cds

clear global lds

format_string = 'Brusselator N:%d  L:%.2f  A:%.0f  B:%.2f  Dx:%.3f  Dy:%.3f\n';
format_args   = {N; L; A; B; Dx; Dy;};
fprintf([format_string '\n'], format_args{:});

initial_continuation_data = init_single_shooting_find_stable_cycle( ...
  'initial_point',             ones(2*N,1), ...
  'time_to_converge_to_cycle', 200, ...
  'subspace_size',             20, ...
  'odefile',                   odefile,  ...
  'ode_parameters',            ode_parameters, ...
  'active_parameter_index',    2, ...
  'lower_bound_period',        1, ...
  'upper_bound_period',        20, ...
  'show_plots',                 true ...
);


fprintf('period: %.4f\n', initial_continuation_data(end-1))

opt = contset( ...
    'MaxNumPoints',                 100, ...
    'InitStepsize',                 1e-1, ...
    'MinStepsize',                  1e-10, ...
    'MaxStepsize',                  5e-2, ...
    'MaxNewtonIters',               3, ...
    'MaxCorrIters',                 6, ...
    'MaxTestIters',                 20, ...
    'VarTolerance',                 1e-6, ...
    'FunTolerance',                 1e-6, ...
    'NewtonPicardBasisTolerance',   1e-6, ...
    'contL_SmoothingAngle',         3, ...
    'CheckClosed',                  1000, ...
    'Multipliers',                  true, ...
    'Backward',                     false, ...
    'Singularities',                true, ...
    'CIS_UsingCIS',                 false, ...
    'NewtonPicard',                 true, ...
    'console_output_level',         4, ...
    'contL_DiagnosticsLevel',       4, ...
    'PicardTolerance',              1e-8);

[s, datafile] = contL(@single_shooting, initial_continuation_data, [], opt); 
