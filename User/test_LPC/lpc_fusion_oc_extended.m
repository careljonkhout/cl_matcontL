% test of limit point bifurcation detection using orhtogonal collocation
clc
N = 25;                     
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
a = -1;
b = -0.3;  
q_inf = -0.72;

ode_parameters = {a ; b; q_inf};

x=loadPoint('/home/carel/Documents/cl_matcontL/User/test_LPC/Data/lpc_fusion_oc_11-Apr-2019_12_17_41.dat');

x(end)

opt = contset();
opt = contset(opt, 'MaxNumPoints',            100);
opt = contset(opt, 'InitStepsize',            0.25);
opt = contset(opt, 'MinStepsize',             1e-6);
opt = contset(opt, 'MaxStepsize',             0.25);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'MaxTestIters',            10);
opt = contset(opt, 'Backward',                true);
opt = contset(opt, 'VarTolerance',            1e-3);
opt = contset(opt, 'FunTolerance',            1e-3);
opt = contset(opt, 'Adapt',                   3);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    3);
opt = contset(opt, 'contL_DiagnosticsLevel',  3);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    1);
opt = contset(opt, 'enable_bpc',              false, ...
                   'enable_nf_pd',            false, ...
                   'enable_nf_lpc',           false, ...
                   'enable_nf_ns',            false);

 

contL(@limitcycleL, x(:,end), [], opt, @plot_T_versus_param);