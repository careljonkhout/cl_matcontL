% continuation of cycles in brusselator using single shooting



function br_ss_cvode



  plotting_enabled = usejava('jvm');



  % N is the number of spatial mesh points in the discretization of the

  % Brusselator model defined in Systems/+Brusselator_1d

  % note that N is fixed by the odefile and the integrator in this case.

  

  N_MESH_POINTS = 100;

  L = 0.6; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;

 

  addpath(fullfile(get_cl_matcontL_path, 'Systems', 'brusselator_1d_cvode'));



  

  if exist('OCTAVE_VERSION', 'builtin') 

    %recompile_for_octave(N_MESH_POINTS);

  else

    recompile_for_matlab(N_MESH_POINTS);

  end

  

  odefile = @odefile_mex;

  integrator = @cvode;

  ode_parameters = {L; A; B; Dx; Dy};

  active_parameter_index = 1;

 

  format_string = 'Brusselator N:%d L:%.2f  A:%.0f  B:%.2f  Dx:%.3f  Dy:%.3f\n';

  format_args   = {N_MESH_POINTS; L; A; B; Dx; Dy;};

  fprintf([format_string '\n'], format_args{:});



  initial_continuation_data = init_single_shooting_find_stable_cycle( ...

    'initial_point',             ones(2 * N_MESH_POINTS, 1), ...

    'time_to_converge_to_cycle', 200, ...

    'subspace_size',             4, ...

    'odefile',                   odefile,  ...

    'time_integration_method',   integrator, ...

    'ode_parameters',            ode_parameters, ...

    'active_parameter_index',    active_parameter_index, ...

    'lower_bound_period',        1, ...

    'upper_bound_period',        20, ...

    'show_plots',                false, ...

    'cvode_verbose',             true ...

  );





  fprintf('period: %.4f\n', initial_continuation_data(end-1))



  opt = contset( ...

      'MaxNumPoints',                 10000, ...

      'InitStepsize',                 5e-2, ...

      'MinStepsize',                  1e-10, ...

      'MaxStepsize',                  5e-3, ...

      'MaxNewtonIters',               3, ...

      'MaxCorrIters',                 10, ...

      'MaxTestIters',                 10, ...

      'FunTolerance',                 10^-8, ...

      'VarTolerance',                 10^-5.5, ...

      'basis_grow_threshold',         0.0001, ...

      'basis_shrink_threshold',       0.000001, ...

      'contL_ParallelComputing',      false, ...

      'num_cores',                    2, ...

      'contL_Testf_FunTolerance',     1e-3, ...

      'contL_Testf_VarTolerance',     1e-3, ...

      'contL_SmoothingAngle',         3, ... % in radians,

      ...                                    % i.e. we disable smoothing by angle

      'CheckClosed',                  1000, ...

      'Multipliers',                  true, ...

      'Backward',                     false, ...

      'Singularities',                true, ...

      'CIS_UsingCIS',                 false, ...

      'NewtonPicard',                 true, ...

      'console_output_level',         2, ...

      'contL_DiagnosticsLevel',       2, ...

      'parameter_sensitivity_by_finite_diff', true, ...

      'PicardTolerance',              1e-13, ...

      'pause',                        false, ...

      'nsteps_before_pause',          1, ...

      'integration_rel_tol',          1e-13, ...

      'integration_abs_tol',          1e-13, ...

      'multipliers_abs_tol',          1e-13, ...

      'multipliers_rel_tol',          1e-13);



  % we open a plot window:

  if plotting_enabled && false

    figure

    % each line segment of the approximation of the curve is plotted separately

    % therefore, the plot has to "hold" the previously plotted line segments

    hold on

    % we set the axes labels

    xlabel('L')

    ylabel('period')

    title_format_string = 'Brusselator N:%d  A:%.0f  B:%.1f  Dx:%.3f  Dy:%.3f';

    title_format_args = {N_MESH_POINTS,A,B,Dx,Dy};

    % we set the title of the plot.

    % the title has two lines

    title({'continuation from stable cycle - single shooting', ...

           sprintf(title_format_string, title_format_args{:})});

  end



  contL(@single_shooting, ...

      initial_continuation_data, [], opt);

end
