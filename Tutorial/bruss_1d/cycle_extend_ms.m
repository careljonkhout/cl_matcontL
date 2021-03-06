% demonstration of how to extend a multiple shooting continuation of a cycle.

% We extend the continuation by loading a .mat file with data from a
% continuation point that was saved during a previous contuation of cycles using
% orthogonal collocation. Try running cycle_integration_ms or
% cycle_from_hopf_ms, and loading one .mat files that one of these scripts saves
% in the "Data" subdirectory instead of loading point_cycle_ms.mat.

% note: the current cl_matcontL implementation (june 2019) of multiple shooting
% with (or without) Newton-Picard is quite slow compared to orthogonal
% collocation for systems of ODEs of up to 350 equations. This file is just
% intended to show how to set up a multiple shooting continuation.

% note: The cycle we will continue here is relatively stable. Therefore, we
% could use single shooting instead of multiple shooting. Single shooting will
% probably always be slightly faster than multiple shooting (regardless of the
% size of the system of ODEs). This file is just intended to show how to set up
% a multiple shooting continuation.

% Note that we define the demo as a function to facilitate testing of the demo
% during the development of cl_matcontL. (Defining a demo as a function prevents
% demo's from being accidentally dependent on some variable in the workspace.)
% It is perfectly fine to run continuations in cl_matcontL using Matlab scripts
% without defining functions.
function cycle_extend_ms

  load(fullfile(get_path(), 'point_cycle_ms.mat'), 'point');


  init_multiple_shooting_extend_curve( ...
    'initial_continuation_state',             point.x, ...
    'n_mesh_intervals',                       point.n_mesh_intervals, ...
    'ode_parameters',                         point.parameter_values, ...
    'active_parameter_index',                 2, ...
    'time_mesh',                              point.time_mesh, ...
    'odefile',                                @brusselator_1d,  ...
    'subspace_size',                          length(point.multipliers));


  opts = contset( ...
    'InitStepsize',           0.05, ...
    'MaxStepsize',            0.05, ...
    'MaxNumPoints',           5, ...
    'contL_DiagnosticsLevel', 4, ...
    'console_output_level',   4, ...
    'NewtonPicard',           true, ...
    'contL_SmoothingAngle',   100, ...
    'newtcorrL_use_max_norm', true, ...
    'Singularities',          true, ...
    'enable_nf_ns',           false, ...
    'enable_nf_lpc',          false, ...
    'enable_nf_pd',           false, ...
    'parameter_sensitivity_by_finite_diff', true, ...
    'singularity_callback',   @plot_singularity_of_cycles);


  % we open a plot window:
  figure
  % each line segment of the approximation of the curve is plotted separately
  % therefore, the plot has to "hold" the previously plotted line segments
  hold on
  % we set the axes labels
  xlabel('L')
  ylabel('period')
  title_format_string = 'Brusselator N:%d  A:%.0f  B:%.1f  Dx:%.3f  Dy:%.3f';
  % we list all the values of the parameters of the system of ODEs in the title,
  % except the active parameter, which is the second parameter in the array of
  % parameters:
  title_format_args = num2cell(point.parameter_values([1,3,4,5,6]));
  % we set the title of the plot.
  % the title has two lines
  title({'extended continuation - multiple shooting', ...
         sprintf(title_format_string, title_format_args{:})});

  contL(@multiple_shooting, point.x, point.v, opts, ...
                      'callback', @plot_T_versus_param);
end


