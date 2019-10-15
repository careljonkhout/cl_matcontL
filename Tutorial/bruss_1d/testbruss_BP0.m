function testbruss_BP0()
% Test script for 1d brusselator (original MATCONT data)

% This example performs continuation of equilibrium curve and locates
% limit points and branch points.

%% Options
opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',             1); 
opt = contset(opt,'contL_DiagnosticsLevel',    3);  
opt = contset(opt,'Backward',                  0);
opt = contset(opt,'InitStepsize',              3);
opt = contset(opt,'MinStepsize',            1e-5);
opt = contset(opt,'MaxStepsize',               3);
opt = contset(opt,'MaxCorrIters',             12);
opt = contset(opt,'MaxNewtonIters',            5);
opt = contset(opt,'FunTolerance',           1e-6);
opt = contset(opt,'VarTolerance',           1e-5);
opt = contset(opt,'contL_SmoothingAngle',   pi/30); 
opt = contset(opt,'Singularities',             1);
opt = contset(opt,'Userfunctions',             0);
opt = contset(opt,'MaxNumPoints',            500);
opt = contset(opt,'MaxTestIters',       10);
%opt = contset(opt,'contL_Testf_VarTolerance', 1e-4); 
%opt = contset(opt,'contL_Testf_FunTolerance', 1e-5); 
opt = contset(opt,'CIS_SparseSolvers',         1);
opt = contset(opt,'CIS_NStableRef',            4);
opt = contset(opt,'CIS_MaxUnstable',           5); %new
opt = contset(opt,'CIS_Ric_Cayley_Shift',      1); 
opt = contset(opt,'contL_EQ_BranchingMethod',  2); 
opt = contset(opt,'Locators',            [1 1 1]);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename',     'testbruss_BP0');

N = 500; L = 0.06; A = 2; B = 4.6; Dx = 0.0016; Dy = 0.008;
p = [N; L; A; B; Dx; Dy]; ap1 = 2;

%% we compute an approximation of an equilibrium x0
x0 = zeros(2*N,1);
i=1:N;
x0(i)   = A + 2*sin(pi*i/(N+1));
x0(N+i) = B/A - 0.5*sin(pi*i/(N+1));

odefile = @brusselator_1d;

%% we initialize the continuation
% init methods initialize the global variable cds, which is a struct contains
% information the continuer needs. 
[x0,v0]       = init_EP_EP_L(odefile, x0, p, ap1);
[s, datafile] = contL(@equilibriumL, x0, v0, opt);

%% Plot results

x = loadPoint(datafile);
N = s(1).data.P0(1);

figure
hold on
title('testbruss\_BP0')
plot(x(2*N+1, :), x(1, :));

for sii = s
    plot(x(2*N+1, sii.index), x(1, sii.index), 'r.');
    text(x(2*N+1, sii.index), x(1, sii.index), sii.label);
end

xlabel 'l'
ylabel 'u(x=0)'