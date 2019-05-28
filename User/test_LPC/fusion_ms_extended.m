clc
clear global
N = 25;                     
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
a = -1;
b = -0.3;
q_inf = -0.72;

ode_parameters = {a ; b; q_inf};



load('/home/carel/Documents/cl_matcontL/User/test_LPC/Data/lpc_fusion_ss_02-Apr-2019_15_27_52.mat')

point = s.data;

init_single_shooting_extend_curve( ...
  'initial_continuation_state',             point.x, ...
  'initial_continuation_tangent_vector',    point.v, ... 
  'ode_parameters',                         num2cell(point.P0), ...
  'active_parameter_index',                 point.ap, ...
  'odefile',                                odefile,  ...
  'subspace_size',                          10 ...
);


opt = contset();
opt = contset(opt, 'MaxNumPoints',            100);
opt = contset(opt, 'InitStepsize',            0.1);
opt = contset(opt, 'MinStepsize',             1e-6);
opt = contset(opt, 'MaxStepsize',             0.1);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'MaxTestIters',            10);
opt = contset(opt, 'Backward',                true);
opt = contset(opt, 'VarTolerance',            1e-6);
opt = contset(opt, 'FunTolerance',            3e-6);
opt = contset(opt, 'Adapt',                   3);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    5);
opt = contset(opt, 'contL_DiagnosticsLevel',  5);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    1);
opt = contset(opt, 'NewtonPicard',            true);
opt = contset(opt, 'integration_rel_tol',              1e-12);
opt = contset(opt, 'integration_abs_tol',              1e-12);
opt = contset(opt, 'multipliers_rel_tol',                1e-13);
opt = contset(opt, 'multipliers_abs_tol',                1e-13);
opt = contset(opt, 'enable_bpc',                       false);
 

contL(@single_shooting, point.x, point.v, opt);