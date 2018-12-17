if ~ exist('contL', 'file')
  fullpath = mfilename('fullpath');
  file_directory = fullpath(1:end-length(mfilename));
  cd([file_directory '../..']);
  init
  cd(file_directory)
end


filename = 'fusion_Orb_LC_10-Dec-2018_17_52_09';
datafile = ['Data/' filename '.dat'];
matrix_file = ['Data/' filename];
load(matrix_file, 's');
load(matrix_file, 'contopts');
%[x, v, h, mult] = loadPoint(datafile);


lpc_point = s(2);
N=25;
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
a = -1;
b = -0.3;
q_inf = lpc_point.data.x(end);
parameter_values = [a,b,q_inf];
ntst=20;
ncol=4;

[x0,v0]= init_LPC_LPC_L(odefile, lpc_point.data.x, lpc_point, [1 3], 20, 4);
opt = contset();
opt = contset(opt, 'InitStepsize',   3);
opt = contset(opt, 'MinStepsize',    1e-8);
opt = contset(opt, 'MaxStepsize',    40);
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
opt = contset(opt, 'Singularities',  true);
opt = contset(opt, 'CIS_UsingCIS',   false);


contL(@limitpointcycle,x0,v0,opt)
