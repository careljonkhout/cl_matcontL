%% Continuation of an equilibrium.
% This example performs continuation of equilibrium curve and locates
% limit points and branch points.

function testbruss_BP0()

%% Options
options = contset( ...
    'contL_DiagnosticsLevel',    3, ...
    'Backward',                  0, ...
    'InitStepsize',              3, ...
    'MinStepsize',               1e-5, ...
    'MaxStepsize',               3, ...
    'MaxCorrIters',              12, ...
    'MaxNewtonIters',            5, ...
    'FunTolerance',              1e-6, ...
    'VarTolerance',              1e-5, ...
    'contL_SmoothingAngle',      pi/30, ... 
    'Singularities',             1, ...
    'Userfunctions',             0, ...
    'MaxNumPoints',              500, ...
    'MaxTestIters',              10, ...
    'CIS_SparseSolvers',         1, ...
    'CIS_NStableRef',            4, ...
    'CIS_MaxUnstable',           5, ...
    'CIS_Ric_Cayley_Shift',      1, ... 
    'contL_EQ_BranchingMethod',  2, ... 
    'Locators',                  [1 1 1]);

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
[singularities, datafile] = contL(@equilibriumL, x0, v0, options);

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