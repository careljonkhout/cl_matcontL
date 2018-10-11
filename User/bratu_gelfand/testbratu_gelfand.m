function testbratu_gelfand()
% Test script for 1d brusselator (original MATCONT data)

% This example locates a user point and saves data for continuation 
% in testbrussL1_U1 with a new parameter

%% Options
opt = contset;
opt = contset(opt,'contL_LogFile',               1);
opt = contset(opt,'contL_DiagnosticsLevel',      3);
opt = contset(opt,'Backward',                    0);
opt = contset(opt,'InitStepsize',             1e-0);
opt = contset(opt,'MinStepsize',              1e-5);
opt = contset(opt,'MaxStepsize',              1e-0);
opt = contset(opt,'MaxCorrIters',               12);
opt = contset(opt,'MaxNewtonIters',              6);
opt = contset(opt,'FunTolerance',             1e-6);
opt = contset(opt,'VarTolerance',             1e-5);
opt = contset(opt,'contL_SmoothingAngle',    pi/30);
opt = contset(opt,'Singularities',               1);
opt = contset(opt,'Userfunctions',               0);
opt = contset(opt,'MaxNumPoints',               30);
opt = contset(opt,'Locators',              [1 1 1]);
opt = contset(opt,'MaxTestIters',                10);
opt = contset(opt,'contL_Testf_VarTolerance', 1e-4);
opt = contset(opt,'contL_Testf_FunTolerance', 1e-5);
opt = contset(opt,'CIS_SparseSolvers',           1);
opt = contset(opt,'CIS_NStableRef',              4);
opt = contset(opt,'CIS_MaxUnstable',             5); %new
opt = contset(opt,'CIS_Ric_Cayley_Shift',       10);
opt = contset(opt,'contL_EQ_BranchingMethod',    2);



%% Continuation


N = 200; LAMBDA = 0;
p = [N; LAMBDA]; ap1 = 2;

[x0,v0]      = init_EP_EP_L(@bratu_gelfand, [], p, ap1);
contL(@equilibriumL,x0,v0,opt);

%{

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
%}