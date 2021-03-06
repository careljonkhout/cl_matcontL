function testfusion_BT0()
% Test script for the fusion model

% This example performs continuation of equilibrium curve and locates a
% limit point.

%% Options
opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',          1);
opt = contset(opt,'contL_DiagnosticsLevel', 4);
opt = contset(opt,'Backward',              0);
opt = contset(opt,'InitStepsize',       0.03);
opt = contset(opt,'MinStepsize',        1e-5);
opt = contset(opt,'MaxStepsize',        0.03);
opt = contset(opt,'Increment',          1e-5);
opt = contset(opt,'FunTolerance',       1e-5);
opt = contset(opt,'VarTolerance',       1e-4);
opt = contset(opt,'contL_SmoothingAngle',  pi/10);
opt = contset(opt,'Singularities',         1);
opt = contset(opt,'Userfunctions',         0);
opt = contset(opt,'MaxNumPoints',         20);
opt = contset(opt,'CIS_SparseSolvers',     1);
opt = contset(opt,'CIS_NStableRef',        6);
opt = contset(opt,'CIS_Ric_Transform',   'cayley');
opt = contset(opt,'CIS_Ric_Cayley_Shift',      4);
opt = contset(opt,'contL_EQ_BranchingMethod',    2);
opt = contset(opt,'Locators',        [1 1 1]);
opt = contset(opt,'contL_Testf_FunTolerance', 1e-4);
opt = contset(opt, 'Filename', 'testfusion_BT0');

%% Parameters

N         =  30;
Gamma_inf = -0.8;
q_inf     = -0.008;
D0        =  1.9;
D1        = -1.1;
D2        =  0;
a         = -0.01;
b         = -0.01; 
zeta1     =  1.1;
mu1       =  0.05;
epsilon   =  0.05; 
ZS        =  0;
gamma1    =  5/3;
lambdan   =  1.25;
lambdaT   =  1.5;
cn        =  1.1;
cT        =  0.9;

ode_parameters = [ N, Gamma_inf, q_inf, D0, D1, D2, a, b, zeta1, mu1, ...
                   epsilon, ZS, gamma1, lambdan, lambdaT, cn, cT];
                 
active_parameter_index = 3;


% initialization
[x0,v0]                   = init_EP_EP_L(@fusion, [], ode_parameters, ...
                                         active_parameter_index); 
% continuation
[singularities, datafile] = contL(@equilibriumL,x0,v0,opt);    

%% Plot results

x = loadPoint(datafile);

figure
hold on
title('testfusion\_BT0')
plot(x(end, :), x(3, :));

for singularity = singularities
  plot(x(end, singularity.index), x(3, singularity.index), 'r.');
  text(x(end, singularity.index), x(3, singularity.index), singularity.label);
end

xlabel 'q_\infty'
ylabel 'Z(x=0)'