%% Continuation of an equilibrium, starting from a point located in testbruss_U0

function testbruss_U1()

%% Options
opt = contset(    'contL_LogFile',               1);
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

%% Load the struct s from the file Data/testbruss_U0.mat
% s contains data about singularities and the first and last points of a
% previous continuation run
path_to_this_script = get_path;
U0_file = fullfile(path_to_this_script, 'Data', 'testbruss_U0.mat');
load(U0_file, 's');

%% We select the third singularity from the data loaded from testbruss_U0.mat
% This is the singularity where the user function specified in the file
% brusselator_1d.m is equal to zero.
ID = 3;
data = s(ID).data;


%% Set A as the active parameter.
% The index of A in the parameter array  ( [N; L; A; B; Dx; Dy] ) is 3. Hence
% the active parameter index must be set to 3 to select A as the active
% parameter.
active_par_idx  = 3;

%% Get the equilibrium from the data that was loaded from testbruss_U0.mat
equilibrium = data.x(1:end-1);

%% Initialize the continuation
% the init functions in cl_matcontL initialize the global variable cds, which is
% a struct contains information the continuer needs.
[x0, v0] = init_EP_EP_L(@brusselator_1d, equilibrium, data.P0, active_par_idx);

%% Continuation
[singularities, datafile] = contL(@equilibriumL, x0, v0, opt);

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