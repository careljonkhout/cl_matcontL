%% Continuation of an equilibrium, with detection of a zero of a "user function"
% The user function is a function supplied by the user in the file
% @brusselator_1d. This example locates the point in the continuation where the
% user function vanishes and saves data for continuation in testbruss_U1 with a
% new parameter.

function testbruss_U0()

%% Options
opt = contset( ...
    'contL_LogFile',               1, ...
    'contL_DiagnosticsLevel',      3, ...
    'Backward',                    1, ...
    'InitStepsize',                1e-0, ...
    'MinStepsize',                 1e-5, ...
    'MaxStepsize',                 1e-0, ...
    'MaxCorrIters',                12, ...
    'MaxNewtonIters',              6, ...
    'FunTolerance',                1e-6, ...
    'VarTolerance',                1e-5, ...
    'contL_SmoothingAngle',        pi/30, ...
    'Singularities',               1, ...
    'Userfunctions',               0, ...
    'MaxNumPoints',                30, ...
    'Locators',                    [1 1 1], ...
    'MaxTestIters',                10, ...
    'contL_Testf_VarTolerance',    1e-4, ...
    'contL_Testf_FunTolerance',    1e-5, ...
    'CIS_SparseSolvers',           1, ...
    'CIS_NStableRef',              4, ...
    'CIS_MaxUnstable',             5, ... 
    'CIS_Ric_Cayley_Shift',        10, ...
    'contL_EQ_BranchingMethod',    2, ...
    'Userfunctions',               1, ...
    'TestPath',                    mfilename('fullpath'));

UserInfo           = [];
UserInfo{1}.name   = 'userf1';
UserInfo{1}.state  = 1;
UserInfo{1}.label  = 'u1';

opt = contset(opt, 'UserFuncInfo', UserInfo);


%% Values of the parameters of the system of ODEs:
N = 128; L = 0.06; A = 2; B = 4.6; Dx = 0.0016; Dy = 0.008;
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
% the init functions in cl_matcontL initialize the global variable cds, which is
% a struct contains information the continuer needs.
[x0,v0] = init_EP_EP_L(@brusselator_1d, equilibrium, p, active_par_idx);

%% Continuation
[singularities, datafile] = contL(@equilibriumL,x0,v0,opt);

%% Plot results


x = loadPoint(datafile);
N = singularities(1).data.P0(1);

figure
hold on
title('testbruss\_U0')
plot(x(N, :), x(N/2, :))

for singularity = singularities
  plot(x(N, singularity.index), x(N/2, singularity.index), 'r.');
  text(x(N, singularity.index), x(N/2, singularity.index), singularity.label);
end

xlabel 'u(x = l)'
ylabel 'u(x = l/2)'