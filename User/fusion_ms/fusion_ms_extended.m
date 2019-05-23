
clc
clear global
N = 25;                     
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
a = -1;
b = -0.3;
q_inf = -0.72;

subdirectory              = sprintf('fusion_ms_%d',N);
dirname                   = fullfile(get_path(), 'Data', subdirectory);
[point_file, point_index] = get_latest_point_file(dirname);
load(point_file, 'point');

ode_parameters = {a ; b; q_inf};


ode_parameters{3} = point.x(end);

init_multiple_shooting_extend_curve( ...
  'initial_continuation_state', point.x, ...
  'time_mesh',                   point.mesh, ...
  'odefile',                     odefile, ...
  'ode_parameters',              ode_parameters, ...
  'active_parameter_index',      3, ...
  'nMeshIntervals',              10, ...
  'subspace_size',               length(point.multipliers) ...
);


opt = contset();
opt = contset(opt, 'MaxNumPoints',            point_index + 30);
opt = contset(opt, 'InitStepsize',            point.h);
opt = contset(opt, 'MinStepsize',             1e-6);
opt = contset(opt, 'MaxStepsize',             0.05);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'MaxTestIters',            10);
opt = contset(opt, 'Backward',                false);
opt = contset(opt, 'VarTolerance',            1e-6);
opt = contset(opt, 'FunTolerance',            1e-6);
opt = contset(opt, 'Adapt',                   3);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    5);
opt = contset(opt, 'contL_DiagnosticsLevel',  5);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    100);
opt = contset(opt, 'NewtonPicard',            true);
opt = contset(opt, 'integration_rel_tol',            1e-9);
opt = contset(opt, 'integration_abs_tol',            1e-9, ...
                   'newtcorrL_use_max_norm',  true, ...
                   'initial_point_index',     point_index, ...
                   'every_point_in_separate_mat_file', true, ...
                   'Filename',                sprintf('fusion_oc_N_%d',N), ...
                   'set_direction',           false);
 

contL(@multiple_shooting, point.x, point.v, opt);