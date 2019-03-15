function testbruss_LP0()
% Test script for 1d brusselator (original MATCONT data)

% This example performs continuation of equilibrium curve and locates 
% a limit point.

%% Options
opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',               1);
opt = contset(opt,'contL_DiagnosticsLevel',      3);
opt = contset(opt,'Backward',                    1);
opt = contset(opt,'InitStepsize',                2);
opt = contset(opt,'MaxStepsize',                 2);
opt = contset(opt,'MinStepsize',              1e-5);
opt = contset(opt,'MaxCorrIters',               12);
opt = contset(opt,'MaxNewtonIters',              5);
opt = contset(opt,'FunTolerance',             1e-6);
opt = contset(opt,'VarTolerance',             1e-5);
opt = contset(opt,'contL_SmoothingAngle',     pi/30);
opt = contset(opt,'Singularities',               1);
opt = contset(opt,'Userfunctions',               0);
opt = contset(opt,'MaxNumPoints',               20);

opt = contset(opt,'Locators',              [1 1 1]);
opt = contset(opt,'contL_Testf_FunTolerance', 1e-6); 
opt = contset(opt,'contL_Testf_VarTolerance', 1e-5);  
opt = contset(opt,'MaxTestIters',               20); 

opt = contset(opt,'CIS_SparseSolvers',           1);
opt = contset(opt,'CIS_NStableRef',              4);
opt = contset(opt,'CIS_MaxUnstable',             5); %new
opt = contset(opt,'CIS_Ric_Cayley_Shift',       10);
opt = contset(opt,'contL_EQ_BranchingMethod',    2);

opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename',     'testbruss_LP0');

%% Continuation 
N = 500; L = 0.06; A = 2; B = 4.6; Dx = 0.0016; Dy = 0.008;
p = [N, L, A, B, Dx, Dy]; ap1 = 2;

[x0,v0]      = init_EP_EP_L(@brusselator_1d, [], p, ap1);
contL(@equilibriumL,x0,v0,opt);

%% Plot results
x = loadPoint('Data\testbruss_LP0.dat');
load('Data\testbruss_LP0.mat')
N = s(1).data.P0(1);
plot(x(2*N+1, :), x(1, :));
hold on
for sii = s
    plot(x(2*N+1, sii.index), x(1, sii.index), 'r.');
    text(x(2*N+1, sii.index), x(1, sii.index), sii.label);
end

xlabel 'l'
ylabel 'u(x=0)'