N = 25;                     
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));

subdirectory              = 'fusion_cycles';
dirname                   = [get_path(), 'Data/', subdirectory];
[point_file, point_index] = get_latest_point_file(dirname);


load(point_file, 'point');

initial_continuation_data = init_collocation_extend_curve( ...
  'continuation_state',           point.x, ...
  'odefile',                      odefile,  ...
  'ode_parameters',               point.parametervalues, ...
  'active_parameter_index',       3, ...
  'time_mesh',                    point.timemesh, ...
  'current_nMeshIntervals',       point.ntst, ...
  'current_nCollocationPoints',   point.ncol ...
);

disp(initial_continuation_data(end-1))
opt = contset();
opt = contset(opt, 'MaxNumPoints',            point_index + 12);
opt = contset(opt, 'InitStepsize',            0.25);
opt = contset(opt, 'MinStepsize',             1e-6);
opt = contset(opt, 'MaxStepsize',             0.25);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'MaxTestIters',            10);
opt = contset(opt, 'VarTolerance',            1e-6);
opt = contset(opt, 'FunTolerance',            1e-6);
opt = contset(opt, 'Adapt',                   3);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    0);
opt = contset(opt, 'contL_DiagnosticsLevel',  0);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    1);
opt = contset(opt,  ...
                   'newtcorrL_use_max_norm',  true, ...
                   'enable_nf_pd',            false, ...
                   'enable_nf_lpc',           true, ...
                   'enable_nf_ns',            false, ...
                   'initial_point_index',     point_index, ...
                   'Filename',                subdirectory, ...
                   'set_direction',           false, ...
                   'singularity_callback',    @plot_singularity_of_cycles);

figure
hold on
title('fusion - extenstion of continuation');
xlabel('q_{inf}')
ylabel('period')
 

contL(@limitcycleL, point.x, point.v, opt,'callback', @plot_T_versus_param);