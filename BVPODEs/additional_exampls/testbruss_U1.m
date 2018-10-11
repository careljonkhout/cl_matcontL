function testbruss_U1()
% Test script for 1d brusselator (original MATCONT data)

% This example performs continuation using active parameter 3 starting from
% a user point located in testbruss_U0
%% Options
opt = contset;
opt = contset(opt,'Cont_LogFile',                1);
opt = contset(opt,'Cont_DiagnosticsLevel',       1);
opt = contset(opt,'Cont_Direction',              0);
opt = contset(opt,'Cont_InitStepsize',        1e-1);
opt = contset(opt,'Cont_MinStepsize',         1e-5);
opt = contset(opt,'Cont_MaxStepsize',         1e-0);
opt = contset(opt,'Cont_MaxCorrIters',          12);
opt = contset(opt,'Cont_MaxNewtonIters',         6);
opt = contset(opt,'Cont_FunTolerance',        1e-6);
opt = contset(opt,'Cont_VarTolerance',        1e-5);
opt = contset(opt,'Cont_SmoothingAngle',     pi/30);
opt = contset(opt,'Cont_Singularities',          1);
opt = contset(opt,'Cont_Userfunctions',          0);
opt = contset(opt,'Cont_MaxNumPoints',          15);
opt = contset(opt,'Loc_UseLocators',       [1 1 1]);
opt = contset(opt,'Loc_Testf_MaxIters',          5);
opt = contset(opt,'Loc_Testf_VarTolerance',   1e-4);
opt = contset(opt,'Loc_Testf_FunTolerance',   1e-5);
opt = contset(opt,'CIS_SparseSolvers',     1);
opt = contset(opt,'CIS_NStableRef',        6);
opt = contset(opt,'RIC_Cayley_Shift',     10);
opt = contset(opt,'EQ_BranchingMethod',    2);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename', 'testbruss_U1');

%% Continuation
load('Data\testbruss_U0.mat')

ID = 3;
data = s(ID).data;
ap1  = 3;
x = data.x(1:end-1);

[x0,v0]      = init_EP_EP_L(@bruss_1d, x, data.P0, ap1);
contL(@equilibriumL,x0,v0,opt);


%% Plot results
x = loadPoint('Data\testbruss_U1.dat');
load('Data\testbruss_U1.mat');
N = s(1).data.P0(1);

hold on
plot(x(N, :), x(N/2, :))
for sii = s
    plot(x(N, sii.index), x(N/2, sii.index), 'r.');
    text(x(N, sii.index), x(N/2, sii.index), sii.label);
end

xlabel 'u(x = l)'
ylabel 'u(x = l/2)'