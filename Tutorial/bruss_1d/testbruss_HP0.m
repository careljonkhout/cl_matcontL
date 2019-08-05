function testbruss_HP0()
% Test script for 1d brusselator (Chien 97 data (!))

% This test function starts from a constant solution of the brussalator and
% detects a hopf point (and possibly branch points that are believed to be
% spurious solutions due to discretizaion.)

%% Options
opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',          1);
opt = contset(opt,'contL_DiagnosticsLevel', 3);
opt = contset(opt,'Backward',               0); %new
opt = contset(opt,'InitStepsize',         0.1);
opt = contset(opt,'MaxStepsize',          0.1);
opt = contset(opt,'MinStepsize',         1e-5);
opt = contset(opt,'Singularities',          1);
opt = contset(opt,'Userfunctions',          0);
opt = contset(opt,'MaxNumPoints',          10);
opt = contset(opt,'Locators',         [1 1 1]);
opt = contset(opt,'CIS_SparseSolvers',      1);
opt = contset(opt,'CIS_MaxUnstable',        5); %new
opt = contset(opt,'CIS_NStableRef',         6);
opt = contset(opt,'CIS_NExtra',             6);
opt = contset(opt,'CIS_Ric_Cayley_Shift',  10);
opt = contset(opt,'contL_Testf_FunTolerance', 1e-10);
opt = contset(opt,'contL_Testf_VarTolerance', 1e-11);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename', 'testbruss_HP0');

%% Start point
N = 200; L = 12; A = 4; B = 17.1; Dx = 1; Dy = 2;
p = [N; L; A; B; Dx; Dy];
ap1 = 4;


%% Continuation
problem_file              = @Brusselator_1d.homogeneous_x0;
[x0,v0]                   = init_EP_EP_L(problem_file, [], p, ap1);
[singularities, datafile] = contL(@equilibriumL,x0,v0,opt);

%% Plot results

x = loadPoint(datafile);
N = singularities(1).data.P0(1);

figure
hold on
title('testbruss\_HP0')
plot(x(2*N+1, :), x(1, :));

for singularity = singularities
  plot(x(2*N+1, singularity.index), x(1, singularity.index), 'r.');
  text(x(2*N+1, singularity.index), x(1, singularity.index), singularity.label);
end

xlabel 'b'
ylabel 'u(x=0)'