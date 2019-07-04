% continuation of cycles in brusselator using single shooting

% note: the current cl_matcontL implementation (june 2019) of single shooting
% with (or without) Newton-Picard is quite slow compared to orthogonal
% collocation for systems of ODEs of up to 350 equations. This file is just
% intended to show how to set up a single shooting continuation.

format long
odefile = @brusselator_1d;
N = 30;
L = 0.5; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
ode_parameters = {N; L; A; B; Dx; Dy};

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
  'show_plots',                true...
);


fprintf('period: %.4f\n', initial_continuation_data(end-1))

opt = contset( ...
    'MaxNumPoints',                 50, ...
    'InitStepsize',                 1e-1, ...
    'MinStepsize',                  1e-10, ...
    'MaxStepsize',                  5e-2, ...
    'MaxNewtonIters',               3, ...
    'MaxCorrIters',                 6, ...
    'MaxTestIters',                 20, ...
    'integration_abs_tol',          1e-9, ...
    'integration_rel_tol',          1e-9, ...
    'FunTolerance',                 1e-5, ...
    'VarTolerance',                 1e-3, ...
    'contL_SmoothingAngle',         3, ... % in radians,
    ...                                    % i.e. we disable smoothing by angle
    'CheckClosed',                  1000, ...
    'Multipliers',                  true, ...
    'Backward',                     false, ...
    'Singularities',                true, ...
    'CIS_UsingCIS',                 false, ...
    'NewtonPicard',                 true, ...
    'console_output_level',         4, ...
    'contL_DiagnosticsLevel',       4, ...
    'PicardTolerance',              1e-8);

% we open a plot window:
figure
% each line segment of the approximation of the curve is plotted separately
% therefore, the plot has to "hold" the previously plotted line segments
hold on
% we set the axes labels
xlabel('L')
ylabel('period')
title_format_string = 'Brusselator N:%d  A:%.0f  B:%.1f  Dx:%.3f  Dy:%.3f';
title_format_args = {N,A,B,Dx,Dy};
% we set the title of the plot.
% the title has two lines
title({'continuation from stable cycle - single shooting', ...
       sprintf(title_format_string, title_format_args{:})});
  
[singularities, datafile] = contL(@single_shooting, ...
    initial_continuation_data, [], opt, 'callback', @plot_T_versus_param);

for singularity = singularities
  plot_singularity_of_cycles(singularity)
end