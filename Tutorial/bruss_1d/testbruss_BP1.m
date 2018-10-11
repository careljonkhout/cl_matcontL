function testbruss_BP1()
% Test script for 1d brusselator (original MATCONT data)

% This example performs continuation of equilibrium curve starting from a
% branch point and continuing along the second branch.

%% Options

opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',          1);
opt = contset(opt,'contL_DiagnosticsLevel', 5);
opt = contset(opt,'Backward',               1);
opt = contset(opt,'InitStepsize',        2e-1);
opt = contset(opt,'MinStepsize',         1e-5);
opt = contset(opt,'MaxStepsize',         1e-0);
opt = contset(opt,'MaxCorrIters',          12);
opt = contset(opt,'MaxNewtonIters',         6);
opt = contset(opt,'FunTolerance',        1e-6);
opt = contset(opt,'VarTolerance',        1e-6);
opt = contset(opt,'contL_SmoothingAngle', pi/30);
opt = contset(opt,'Singularities',          1);
opt = contset(opt,'Userfunctions',          0);
opt = contset(opt,'MaxNumPoints',         500);
opt = contset(opt,'CIS_SparseSolvers',      1);
opt = contset(opt,'CIS_NStableRef',         6);
opt = contset(opt,'CIS_MaxUnstable',        5); %new
opt = contset(opt,'CIS_Ric_Cayley_Shift',  10);
opt = contset(opt,'contL_EQ_BranchingMethod',2);
opt = contset(opt,'Locators',        [1 1 1]);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename', 'testbruss_BP1');

%% Continuation
load('Data\testbruss_BP0.mat')

ID = 3;
if strcmp(s(ID).label,'BP')
    data = s(ID).data;
    [x0,v0]      = init_BP_EP_L(@bruss_1d, [], [], [], data);
    contL(@equilibriumL,x0,v0,opt);
else
    debug('No Branch Vector\n');
end

%% Plot results

x = loadPoint('Data\testbruss_BP1.dat');
load('Data\testbruss_BP1.mat');
hold on
N = s(1).data.P0(1);
plot(x(2*N+1, :), x(1, :));
for sii = s
    plot(x(2*N+1, sii.index), x(1, sii.index), 'r.');
    text(x(2*N+1, sii.index), x(1, sii.index), sii.label);
end

xlabel 'l'
ylabel 'u(x=0)'