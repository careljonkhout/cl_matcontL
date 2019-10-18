%% Continuation of equilibrium, with detection of a Hopf point
% This example starts from a constant solution of the Brussalator and detects a
% hopf point (and possibly branch points that are believed to be spurious
% solutions due to discretizaion.)

function testbruss_HP0()

%% Options
opt = contset(    'contL_LogFile',          1);
opt = contset(opt,'contL_DiagnosticsLevel', 3);
opt = contset(opt,'Backward',               0);
opt = contset(opt,'InitStepsize',         0.1);
opt = contset(opt,'MaxStepsize',          0.1);
opt = contset(opt,'MinStepsize',         1e-5);
opt = contset(opt,'Singularities',          1);
opt = contset(opt,'Userfunctions',          0);
opt = contset(opt,'MaxNumPoints',          10);
opt = contset(opt,'Locators',         [1 1 1]);
opt = contset(opt,'CIS_SparseSolvers',      1);
opt = contset(opt,'CIS_MaxUnstable',        5);
opt = contset(opt,'CIS_NStableRef',         6);
opt = contset(opt,'CIS_NExtra',             6);
opt = contset(opt,'CIS_Ric_Cayley_Shift',  10);
opt = contset(opt,'contL_Testf_FunTolerance', 1e-10);
opt = contset(opt,'contL_Testf_VarTolerance', 1e-11);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename', 'testbruss_HP0');

%% Values of the parameters of the system of ODEs:
N = 200; L = 12; A = 4; B = 17.1; Dx = 1; Dy = 2;
p = [N; L; A; B; Dx; Dy];

%% Use Dx as the active parameter.
% The index of Dx in the parameter array is 4. Hence the active parameter index
% must be set to 4 to select Dx as the active parameter.
active_par_idx = 4;

%% Compute the initial point "equilibrium".
% This equilibrium represents a spatially homogeneous equilibrium of the PDE:
equilibrium          = zeros(2*N,1);
equilibrium(1:N)     = A;
equilibrium(N+1:2*N) = B/A;

%% Initialize the continuation
% the init functions in cl_matcontL initialize the global variable cds, which is
% a struct contains information the continuer needs.
[x0,v0] = init_EP_EP_L(@brusselator_1d, equilibrium, p, active_par_idx);

%% Continuation
[singularities, datafile] = contL(@equilibriumL, x0, v0, opt);

%% Plot results

x = loadPoint(datafile);

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