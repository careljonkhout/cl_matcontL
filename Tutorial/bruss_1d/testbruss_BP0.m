%% Continuation of an equilibrium.
% This example performs continuation of equilibrium curve and locates
% limit points and branch points.

function testbruss_BP0()

%% Options
opt = contset(    'contL_LogFile',             1); 
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

%% Values of the parameters of the system of ODEs:
N = 500; L = 0.06; A = 2; B = 4.6; Dx = 0.0016; Dy = 0.008;
p = [N; L; A; B; Dx; Dy];

%% Use the L as the active parameter.
% The index of L in the parameter array is 2. Hence the active parameter index
% must be set to 2 to select L as the active parameter.
active_par_idx = 2;

%% Compute an approximation of an equilibrium
% "equilibrium" represents a spatially heterogenous equilibrium of the PDE
equilibrium      = zeros(2*N,1);
i                = 1:N;
equilibrium(  i) = A   + 2   * sin( pi * i / (N+1) );
equilibrium(N+i) = B/A - 0.5 * sin( pi * i / (N+1) );

%% Initialize the continuation
% init methods initialize the global variable cds, which is a struct contains
% information the continuer needs. 
[x0,v0]       = init_EP_EP_L(@brusselator_1d, equilibrium, p, active_par_idx);
%% Continuation
[singularities, datafile] = contL(@equilibriumL, x0, v0, opt);

%% Plot results

x = loadPoint(datafile);

figure
hold on
title('testbruss\_BP0')
plot(x(2*N+1, :), x(1, :));

for sii = singularities
    plot(x(2*N+1, sii.index), x(1, sii.index), 'r.');
    text(x(2*N+1, sii.index), x(1, sii.index), sii.label);
end

xlabel 'l'
ylabel 'u(x=0)'