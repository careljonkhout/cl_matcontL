function testbruss_U0()
% Test script for 1d brusselator (original MATCONT data)

% This example locates a user point and saves data for continuation 
% in testbrussL1_U1 with a new parameter

%% Options
opt = contset;
opt = contset(opt,'Cont_LogFile',                1);
opt = contset(opt,'Cont_DiagnosticsLevel',       3);
opt = contset(opt,'Cont_Direction',              1);
opt = contset(opt,'Cont_InitStepsize',        1e-0);
opt = contset(opt,'Cont_MinStepsize',         1e-5);
opt = contset(opt,'Cont_MaxStepsize',         1e-0);
opt = contset(opt,'Cont_MaxCorrIters',          12);
opt = contset(opt,'Cont_MaxNewtonIters',         6);
opt = contset(opt,'Cont_FunTolerance',        1e-6);
opt = contset(opt,'Cont_VarTolerance',        1e-5);
opt = contset(opt,'Cont_SmoothingAngle',     pi/30);
opt = contset(opt,'Cont_Singularities',          1);
opt = contset(opt,'Cont_Userfunctions',          0);
opt = contset(opt,'Cont_MaxNumPoints',          30);
opt = contset(opt,'Loc_UseLocators',       [1 1 1]);
opt = contset(opt,'Loc_Testf_MaxIters',          5);
opt = contset(opt,'Loc_Testf_VarTolerance',   1e-4);
opt = contset(opt,'Loc_Testf_FunTolerance',   1e-5);
opt = contset(opt,'CIS_SparseSolvers',           1);
opt = contset(opt,'CIS_NStableRef',              4);
opt = contset(opt,'RIC_Cayley_Shift',           10);
opt = contset(opt,'EQ_BranchingMethod',          2);
opt = contset(opt,'Cont_Userfunctions',          1);

opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename',      'testbruss_U0');

UserInfo           = [];
UserInfo{1}.name   = 'userf1';
UserInfo{1}.state  = 1;
UserInfo{1}.label  = 'u1';
opt = contset(opt,'Cont_UserFuncInfo',    UserInfo);

%% Continuation
N = 128; L = 0.06; A = 2; B = 4.6; Dx = 0.0016; Dy = 0.008;
p = [N; L; A; B; Dx; Dy]; ap1 = 2;

[x0,v0]      = init_EP_EP_L(@bruss_1d, [], p, ap1);
contL(@equilibriumL,x0,v0,opt);

%% Plot results
x = loadPoint('Data\testbruss_U0.dat');
load('Data\testbruss_U0.mat');

plot(x(N, :), x(N/2, :))
hold on
for sii = s
    plot(x(N, sii.index), x(N/2, sii.index), 'r.');
    text(x(N, sii.index), x(N/2, sii.index), sii.label);
end

xlabel 'u(x = l)'
ylabel 'u(x = l/2)'