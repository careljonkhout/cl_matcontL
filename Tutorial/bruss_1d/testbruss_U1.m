function testbruss_U1()
% Test script for 1d brusselator (original MATCONT data)

% This example performs continuation using active parameter 3 starting from
% a user point located in testbruss_U0
%% Options
opt = contset;
opt = contset(opt,'contL_LogFile',               1);
opt = contset(opt,'contL_DiagnosticsLevel',      1);
opt = contset(opt,'Backward',                    0);
opt = contset(opt,'InitStepsize',             1e-1);
opt = contset(opt,'MinStepsize',              1e-5);
opt = contset(opt,'MaxStepsize',              1e-0);
opt = contset(opt,'MaxCorrIters',               12);
opt = contset(opt,'MaxNewtonIters',              6);
opt = contset(opt,'FunTolerance',             1e-6);
opt = contset(opt,'VarTolerance',             1e-5);
opt = contset(opt,'contL_SmoothingAngle',    pi/30);
opt = contset(opt,'Singularities',               1);
opt = contset(opt,'Userfunctions',               0);
opt = contset(opt,'MaxNumPoints',               15);
opt = contset(opt,'Locators',              [1 1 1]);
opt = contset(opt,'MaxTestIters',                5);
opt = contset(opt,'contL_Testf_VarTolerance', 1e-4);
opt = contset(opt,'contL_Testf_FunTolerance', 1e-5);
opt = contset(opt,'CIS_SparseSolvers',           1);
opt = contset(opt,'CIS_NStableRef',              6);
opt = contset(opt,'CIS_MaxUnstable',             5); %new
opt = contset(opt,'CIS_Ric_Cayley_Shift',       10);
opt = contset(opt,'contL_EQ_BranchingMethod',    2);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename', 'testbruss_U1');

%% Continuation
path_to_this_script = get_path;
U0_file = [path_to_this_script, 'Data/testbruss_U0.mat'];
load(U0_file, 's');


ID = 3;
data = s(ID).data;
ap1  = 3;
x = data.x(1:end-1);

[x0,v0]      = init_EP_EP_L(@brusselator_1d, x, data.P0, ap1);
[singularities, datafile] = contL(@equilibriumL,x0,v0,opt);

x = loadPoint(datafile);
N = singularities(1).data.P0(1);

figure
hold on
title('testbruss\_U1')
plot(x(N, :), x(N/2, :))

for singularity = singularities
  plot(x(N, singularity.index), x(N/2, singularity.index), 'r.');
  text(x(N, singularity.index), x(N/2, singularity.index), singularity.label);
end

xlabel 'u(x = l)'
ylabel 'u(x = l/2)'