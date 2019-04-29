% test of limit point bifurcation detection using orhtogonal collocation
clc
clear global
N = 75;                     
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
a = -1;
b = -0.3;  
q_inf = -0.72;

ode_parameters = {a ; b; q_inf};


initial_continuation_data = init_collocation_extend_curve( ...
  'continuation_state',           point.x, ...
  'continuation_tangent',         point.v, ... 
  'odefile',                      odefile,  ...
  'ode_parameters',               point.parametervalues, ...
  'active_parameter_index',       3, ...
  'time_mesh',                    point.timemesh, ...
  'current_nMeshIntervals',       40, ...
  'current_nCollocationPoints',   4 ...
);

disp(initial_continuation_data(end-1))
opt = contset();
opt = contset(opt, 'MaxNumPoints',            20);
opt = contset(opt, 'InitStepsize',            0.25);
opt = contset(opt, 'MinStepsize',             1e-6);
opt = contset(opt, 'MaxStepsize',             0.25);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'MaxTestIters',            10);
opt = contset(opt, 'Backward',                false);
opt = contset(opt, 'VarTolerance',            1e-6);
opt = contset(opt, 'FunTolerance',            1e-6);
opt = contset(opt, 'Adapt',                   3);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    3);
opt = contset(opt, 'contL_DiagnosticsLevel',  3);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    1);
opt = contset(opt,  ...
                   'enable_nf_pd',            false, ...
                   'enable_nf_lpc',           true, ...
                   'enable_nf_ns',            false);

 

contL(@limitcycleL, point.x, point.v, opt);