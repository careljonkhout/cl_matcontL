filename = 'fusion_Orb_LC_10-Dec-2018_17_52_09';
matrix_file = fullfile(get_path(), 'Data', 'old_data', filename);
load(matrix_file, 's');
load(matrix_file, 'contopts');
%[x, v, h, mult] = loadPoint(datafile);


lpc_point = s(2);
N=25;
%odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
odefile = str2func(sprintf('fusion_N_%d_hess', N));
a = -1;
b = -0.3;
q_inf = lpc_point.data.x(end);
parameter_values = [a,b,q_inf];


[x0,v0]= init_LPC_LPC_L(odefile, lpc_point.data.x, lpc_point, [1 3], 80, 4);
opt = contset();
opt = contset(opt, 'InitStepsize',   1e-1);
opt = contset(opt, 'MinStepsize',    1e-8);
opt = contset(opt, 'MaxStepsize',    1e-1);
opt = contset(opt, 'MaxNewtonIters', 6);
opt = contset(opt, 'MaxCorrIters',   18);
opt = contset(opt, 'MaxTestIters',   10);
opt = contset(opt, 'VarTolerance',   10^(-2.5));
opt = contset(opt, 'FunTolerance',   5e-5);
opt = contset(opt, 'Adapt',          3); 
opt = contset(opt, 'MaxNumPoints',   100000, ...
      'console_output_level',             5, ...
    'contL_DiagnosticsLevel',           5, ...
                    'MoorePenrose',  false);
opt = contset(opt, 'CheckClosed',    50);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       false);
opt = contset(opt, 'Singularities',  false);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'CIS_UsingCIS',   false);
opt = contset(opt, 'newtcorrL_use_max_norm', true);
opt = contset(opt, 'contL_SmoothingAngle', pi/2/10);

contL(@limitpointcycle,x0,v0,opt)
