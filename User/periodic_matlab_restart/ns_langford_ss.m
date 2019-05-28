clc
clear global

dirname  = fullfile(get_path(),'Data',mfilename);
filename = get_latest_mat_file(dirname);
initial_point_index = sscanf(filename,'point_%d');
load(fullfile(dirname,filename));

init_single_shooting_extend_curve( ...
  'initial_continuation_state',             point.x, ...
  'ode_parameters',                         num2cell([3,0.25,0.2,1.9]), ...
  'active_parameter_index',                 4, ...
  'odefile',                                @Langford,  ...
  'subspace_size',                          10 ...
);


opt = contset();
opt = contset(opt, 'MaxNumPoints',             initial_point_index + 3);
opt = contset(opt, 'InitStepsize',             4e-2);
opt = contset(opt, 'MinStepsize',              1e-10);
opt = contset(opt, 'MaxStepsize',              4e-2);
opt = contset(opt, 'MaxNewtonIters',           8);
opt = contset(opt, 'MaxCorrIters',             10);
opt = contset(opt, 'VarTolerance',             1e-5);
opt = contset(opt, 'FunTolerance',             1e-9);
opt = contset(opt, 'MaxTestIters',             40);
opt = contset(opt, 'contL_Testf_VarTolerance', 1e-6);
opt = contset(opt, 'contL_Testf_FunTolerance', 1e-6);
opt = contset(opt, 'Adapt',                    1000*1000*1000);
opt = contset(opt, 'Multipliers',              true);
opt = contset(opt, 'Singularities',            true);
opt = contset(opt, 'console_output_level',     1);
opt = contset(opt, 'contL_DiagnosticsLevel',   0);
opt = contset(opt, 'MoorePenrose',             false);
opt = contset(opt, 'contL_SmoothingAngle',     pi, ...
                   'pause',                    true, ...
                   'nsteps_before_pause',      20, ...
                   'Filename',                 mfilename, ... 
                   'initial_point_index',               initial_point_index);

global contopts
contopts = opt;


contL(@single_shooting, point.x, point.v, opt);