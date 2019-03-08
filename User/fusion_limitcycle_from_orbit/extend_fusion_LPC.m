run_init_if_needed

filename = 'fusion_LPC_17-Feb-2019_9_40_03';
datafile = ['Data/' filename '.dat'];
matrix_file = ['Data/' filename];
load(matrix_file, 's');
load(matrix_file, 'contopts');
%[x, v, h, mult] = loadPoint(datafile);


last_point = s(end);
N=25;
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
a = -1;
b = -0.3;
q_inf = last_point.data.x(end);
parameter_values = [a,b,q_inf];
ntst=80;
ncol=4;

[x0,v0]= init_LPC_LPC_L(odefile, last_point.data.x, last_point, [1 3], 20, 4);
opt = contset();
opt = contset(opt, 'InitStepsize',   1);
opt = contset(opt, 'MinStepsize',    1e-8);
opt = contset(opt, 'MaxStepsize',    1);
opt = contset(opt, 'MaxNewtonIters', 6);
opt = contset(opt, 'MaxCorrIters',   18);
opt = contset(opt, 'MaxTestIters',   10);
opt = contset(opt, 'VarTolerance',   1e-5);
opt = contset(opt, 'FunTolerance',   1e-5);
opt = contset(opt, 'Adapt',          3); 
opt = contset(opt, 'MaxNumPoints',   1000);
opt = contset(opt, 'CheckClosed',    50);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       false);
opt = contset(opt, 'Singularities',  false);
opt = contset(opt, 'CIS_UsingCIS',   false);
opt = contset(opt, 'newtcorrL_use_max_norm', true);
% disable smoothing by angle:
opt = contset(opt, 'contL_SmoothingAngle', pi/2);

contL(@limitpointcycle,x0,v0,opt)
