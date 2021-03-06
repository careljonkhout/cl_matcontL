% continuation of cycles in brusselator using orthogonal collocation from a
% stable cycle found by time integration

% Note that we define the demo as a function to facilitate testing of the demo
% during the development of cl_matcontL. (Defining a demo as a function prevents
% demo's from being accidentally dependent on some variable in the workspace.)
% It is perfectly fine to run continuations in cl_matcontL using Matlab scripts
% without defining functions.
function cycle_integration_oc
  N = 30;
  L = 0.52; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
  ode_parameters = {N; L; A; B; Dx; Dy};

  format_string = 'Brusselator N:%d  L:%.2f  A:%.0f  B:%.2f  Dx:%.3f  Dy:%.3f\n';
  format_args   = {N; L; A; B; Dx; Dy;};
  fprintf([format_string '\n'], format_args{:});

  initial_continuation_data = init_collocation_find_stable_cycle( ...
    'initial_point',             ones(2*N,1), ...
    'time_to_converge_to_cycle', 200, ...
    'odefile',                   @brusselator_1d,  ...
    'ode_parameters',            ode_parameters, ...
    'n_mesh_intervals',            20, ...
    'active_parameter_index',    2, ...
    'lower_bound_period',        1, ...
    'upper_bound_period',        4, ...
    'poincare_tolerance',        0.1, ...
    'show_plots',                true ...
  );


  fprintf('period: %.4f\n', initial_continuation_data(end-1))

  opt = contset( ...
      'MaxNumPoints',                 70, ...
      'InitStepsize',                 1e-1, ...
      'MinStepsize',                  1e-10, ...
      'MaxStepsize',                  5e-2, ...
      'MaxNewtonIters',               3, ...
      'MaxCorrIters',                 6, ...
      'MaxTestIters',                 20, ...
      'VarTolerance',                 1e-6, ...
      'FunTolerance',                 1e-6, ...
      'contL_SmoothingAngle',         3, ...
      'CheckClosed',                  1000, ...
      'Multipliers',                  true, ...
      'Backward',                     false, ...
      'Singularities',                true, ...
      'CIS_UsingCIS',                 false, ...
      'console_output_level',         0, ...
      'contL_DiagnosticsLevel',       0, ...
      'newtcorrL_use_max_norm',       true, ...
      'PicardTolerance',              1e-8, ...
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
  title({'continuation from stable cycle - orthogonal collocation', ...
         sprintf(title_format_string, title_format_args{:})});

  contL(@limitcycleL, ...
    initial_continuation_data, [], opt, 'callback', @plot_T_versus_param); 
end
