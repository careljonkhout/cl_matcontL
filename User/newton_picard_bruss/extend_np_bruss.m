load('/home/carel/Documents/cl_matcontL/User/newton_picard_bruss/Data/newton_picard_bruss_19-Apr-2019_15_24_52/point 49.mat')



N = (length(point.x)-2)/2;
L = 1.1; A = 2; B = 5.45; Dx = 0.008; Dy = 0.004;
ode_parameters = {N; L; A; B; Dx; Dy};%parameters = {L; A; B; Dx; Dy};
ode_parameters{2} = point.x(end);
init_single_shooting_extend_curve(...
 'initial_continuation_state', point.x, ...
 'initial_continuation_tangent_vector', point.v, ...
 'odefile',  @brusselator_1d, ...
 'ode_parameters', ode_parameters, ...
 'active_parameter_index', 2, ...
 'subspace_size', length(point.multipliers));


opt = contset();
opt = contset(opt, 'MaxNumPoints',   10000);
opt = contset(opt, 'InitStepsize',   1e-1);
opt = contset(opt, 'MinStepsize',    1e-10);
opt = contset(opt, 'MaxStepsize',    5e-2);
opt = contset(opt, 'MaxNewtonIters', 3);
opt = contset(opt, 'MaxCorrIters',   6);
opt = contset(opt, 'MaxTestIters',   20);
opt = contset(opt, 'VarTolerance',   1e-6);
opt = contset(opt, 'FunTolerance',   1e-6);
opt = contset(opt, 'NewtonPicardBasisTolerance',   1e-6);
opt = contset(opt, 'contL_SmoothingAngle',   3);
% we don't want to adapt
% since it is not implemented
opt = contset(opt, 'Adapt',          1000*1000*1000);
opt = contset(opt, 'CheckClosed',    1000);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       false);
opt = contset(opt, 'Singularities',  true);
opt = contset(opt, 'CIS_UsingCIS',   false);
opt = contset(opt, 'NewtonPicard',   true);
opt = contset(opt, 'console_output_level',   4);
opt = contset(opt, 'contL_DiagnosticsLevel', 4);
opt = contset(opt, 'every_point_in_separate_mat_file', true);
opt = contset(opt, 'PicardTolerance', 1e-8);

s = contL(@single_shooting,point.x,[],opt);