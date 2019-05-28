% test of Neimark Sacker bifurcation detection using orhtogonal collocation
clc
clear global
N = 100;                     
odefile = @nonadiabatic_tubular_reactor;
D     = 0.165;
P_em  = 5;
P_eh  = 5;
BETA  = 2.35;
phi_0 = 1;
GAMMA = 25;
B     = 0.5;

ode_parameters = {D; P_em; P_eh; BETA; phi_0; GAMMA; B};


opt = contset();
opt = contset(opt, 'CIS_UsingCIS',            false);
opt = contset(opt, 'MaxNumPoints',            100);
opt = contset(opt, 'InitStepsize',            1);
opt = contset(opt, 'MinStepsize',             1e-6);
opt = contset(opt, 'MaxStepsize',             1);
opt = contset(opt, 'MaxNewtonIters',          8);
opt = contset(opt, 'MaxCorrIters',            10);
opt = contset(opt, 'MaxTestIters',            10);
opt = contset(opt, 'Backward',                false);
opt = contset(opt, 'VarTolerance',            1e-3);
opt = contset(opt, 'FunTolerance',            1e-3);
opt = contset(opt, 'Adapt',                   3);
opt = contset(opt, 'Multipliers',             true);
opt = contset(opt, 'Singularities',           true);
opt = contset(opt, 'console_output_level',    3);
opt = contset(opt, 'contL_DiagnosticsLevel',  3);
opt = contset(opt, 'MoorePenrose',            false);
opt = contset(opt, 'contL_SmoothingAngle',    1);
opt = contset(opt, 'enable_bpc',              false);

 
[x,v] = init_EP_EP_L(odefile,ones(2*N,1),cell2mat(ode_parameters),1);
contL(@equilibriumL, x,v, opt);