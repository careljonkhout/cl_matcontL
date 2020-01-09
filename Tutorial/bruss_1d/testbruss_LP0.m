%% Continuation of equilibrium, with detection of a limit point
% This example performs continuation of equilibrium curve and locates 
% a limit point in the 1d brusselator (Chien 97 data)

function testbruss_LP0()


%% Options
opt = contset( ...
    'contL_DiagnosticsLevel',      3, ...
    'Backward',                    1, ...
    'InitStepsize',                2, ...
    'MaxStepsize',                 2, ...
    'MinStepsize',                 1e-5, ...
    'MaxCorrIters',                12, ...
    'MaxNewtonIters',              5, ...
    'FunTolerance',                1e-6, ...
    'VarTolerance',                1e-5, ...
    'contL_SmoothingAngle',        pi/30, ...
    'Singularities',               1, ...
    'Userfunctions',               0, ...
    'MaxNumPoints',                20, ...
    'Locators',                    [1 1 1], ...
    'contL_Testf_FunTolerance',    1e-6, ... 
    'contL_Testf_VarTolerance',    1e-5, ...  
    'MaxTestIters',                20, ... 
    'CIS_SparseSolvers',           1, ...
    'CIS_NStableRef',              4, ...
    'CIS_MaxUnstable',             5, ... %new
    'CIS_Ric_Cayley_Shift',        10, ...
    'contL_EQ_BranchingMethod',    2, ...
    'TestPath',                    mfilename('fullpath'));


%% Values of the parameters of the system of ODEs:
N = 500; L = 0.06; A = 2; B = 4.6; Dx = 0.0016; Dy = 0.008;
p = [N, L, A, B, Dx, Dy];

%% Use L as the active parameter.
% The index of L in the parameter array is 2. Hence the active parameter index
% must be set to 2 to select L as the active parameter.
active_par_idx = 2;


%% Compute an approximation of an equilibrium x0
% x0 represents a spatially heterogenous equilibrium of the PDE
equilibrium      = zeros(2*N,1);
i                = 1:N;
equilibrium(  i) = A   + 2   * sin( pi * i / (N+1) );
equilibrium(N+i) = B/A - 0.5 * sin( pi * i / (N+1) );

%% Initialize the continuation
% the init functions in cl_matcontL initialize the global variable cds, which is
% a struct contains information the continuer needs.
[x0, v0] = init_EP_EP_L(@brusselator_1d, equilibrium, p, active_par_idx);

%% Continuation 
[singularities, datafile] = contL(@equilibriumL, x0, v0, opt);

x = loadPoint(datafile);
N = singularities(1).data.P0(1);

figure
hold on
title('testbruss\_LP0')
plot(x(2*N+1, :), x(1, :));

for singularity = singularities
  plot(x(2*N+1, singularity.index), x(1, singularity.index), 'r.');
  text(x(2*N+1, singularity.index), x(1, singularity.index), singularity.label);
end

xlabel 'l'
ylabel 'u(x=0)'