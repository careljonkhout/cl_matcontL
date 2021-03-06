function testbruss_BP0()
% Test script for 1d brusselator (original MATCONT data)

% This example performs continuation of equilibrium curve and locates
% limit points and branch points.

%% Options
opt = contset(); %Clear previous options
opt = contset(opt,'Cont_LogFile',              1);
opt = contset(opt,'Cont_DiagnosticsLevel',     3);
opt = contset(opt,'Cont_Direction',            0);
opt = contset(opt,'Cont_InitStepsize',         3);
opt = contset(opt,'Cont_MinStepsize',       1e-5);
opt = contset(opt,'Cont_MaxStepsize',          3);
opt = contset(opt,'Cont_MaxCorrIters',        12);
opt = contset(opt,'Cont_MaxNewtonIters',       5);
opt = contset(opt,'Cont_FunTolerance',      1e-6);
opt = contset(opt,'Cont_VarTolerance',      1e-5);
opt = contset(opt,'Cont_SmoothingAngle',   pi/30);
opt = contset(opt,'Cont_Singularities',        1);
opt = contset(opt,'Cont_Userfunctions',        0);
opt = contset(opt,'Cont_MaxNumPoints',       500);
opt = contset(opt,'Loc_Testf_MaxIters',       10);
opt = contset(opt,'Loc_Testf_VarTolerance', 1e-4);
opt = contset(opt,'Loc_Testf_FunTolerance', 1e-5);
opt = contset(opt,'CIS_SparseSolvers',         1);
opt = contset(opt,'CIS_NStableRef',            4);
opt = contset(opt,'RIC_Cayley_Shift',         1);
opt = contset(opt,'EQ_BranchingMethod',        2);
opt = contset(opt,'Loc_UseLocators',     [1 1 1]);

opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename',     'testbruss_BP0');

%% Continuation
N = 500; L = 0.06; A = 2; B = 4.6; Dx = 0.0016; Dy = 0.008;
p = [N; L; A; B; Dx; Dy]; ap1 = 2;

[x0,v0]      = init_EP_EP_L(@bruss_1d, [], p, ap1);
contL(@equilibriumL,x0,v0,opt);

%% Plot results

x = loadPoint('Data\testbruss_BP0.dat');
load('Data\testbruss_BP0.mat');
N = s(1).data.P0(1);
plot(x(2*N+1, :), x(1, :));
hold on
for sii = s
    plot(x(2*N+1, sii.index), x(1, sii.index), 'r.');
    text(x(2*N+1, sii.index), x(1, sii.index), sii.label);
end

xlabel 'l'
ylabel 'u(x=0)'