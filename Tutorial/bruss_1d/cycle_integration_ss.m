% continuation of cycles in brusselator using single shooting

% note: the current cl_matcontL implementation (june 2019) of single shooting
% with (or without) Newton-Picard is quite slow compared to orthogonal
% collocation for systems of ODEs of up to 350 equations. This file is just
% intended to show how to set up a single shooting continuation.

% Note that we define the demo as a function to facilitate testing of the demo
% during the development of cl_matcontL. (Defining a demo as a function prevents
% demo's from being accidentally dependent on some variable in the workspace.)
% It is perfectly fine to run continuations in cl_matcontL using Matlab scripts
% without defining functions.
function cycle_integration_ss

  % N is the number of spatial mesh points in the discretization of the
  % Brusselator model defined in Systems/+Brusselator_1d
  N = 30;
  L = 0.5; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
  ode_parameters = {N; L; A; B; Dx; Dy};

  format_string = 'Brusselator N:%d L:%.2f A:%.0f B:%.2f Dx:%.3f Dy:%.3f';
  format_args   = {N; L; A; B; Dx; Dy;};
  fprintf([format_string '\n'], format_args{:});

  initial_continuation_data = init_single_shooting_find_stable_cycle( ...
    'initial_point',             ones(2*N,1), ...
    'time_to_converge_to_cycle', 200, ...
    'subspace_size',             20, ...
    'odefile',                   @brusselator_1d,  ...
    'ode_parameters',            ode_parameters, ...
    'active_parameter_index',    2, ...
    'lower_bound_period',        1, ...
    'upper_bound_period',        20, ...
    'show_plots',                true...
  );


  fprintf('period: %.4f\n', initial_continuation_data(end-1))

  opt = contset( ...
      'MaxNumPoints',                 23, ...
      'InitStepsize',                 5e-2, ...
      'MinStepsize',                  1e-10, ...
      'MaxStepsize',                  5e-2, ...
      'MaxNewtonIters',               3, ...
      'MaxCorrIters',                 6, ...
      'MaxTestIters',                 20, ...
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
      'PicardTolerance',              1e-8, ...
      'parameter_sensitivity_by_finite_diff', true, ...
      'singularity_callback',         @plot_singularity_of_cycles);

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

  contL(@single_shooting, ...
      initial_continuation_data, [], opt, 'callback', @plot_T_versus_param);
end
