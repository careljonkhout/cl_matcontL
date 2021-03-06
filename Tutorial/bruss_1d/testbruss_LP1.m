%% Continuation of a limit point
% This example performs a continuation of a limit point starting from a point
% located in testbruss_LP0

function testbruss_LP1()

%% Options
opt = contset( ...
    'contL_LogFile',              1, ...
    'contL_DiagnosticsLevel',     5, ...
    'Backward',                   0, ...
    'InitStepsize',               3, ...
    'MaxStepsize',                3, ...
    'MinStepsize',                1e-8, ...
    'MaxCorrIters',               12, ...
    'MaxNewtonIters',             5, ...
    'FunTolerance',               1e-5, ...
    'VarTolerance',               1e-5, ...
    'contL_SmoothingAngle',       pi/30, ...
    'Singularities',              1, ...
    'Userfunctions',              0, ...
    'MaxNumPoints',               150, ...
    'contL_Testf_FunTolerance',   1e-4, ...
    'contL_Testf_VarTolerance',   1e-4, ...
    'MaxTestIters',               16, ...
    'CIS_SparseSolvers',          1, ...
    'CIS_NStableRef',             9, ...
    'CIS_MaxUnstable',            5, ... 
    'CIS_NExtra',                 4, ...
    'CIS_Ric_Cayley_Shift',       10);

%% Load the struct s from the file Data/testbruss_LP0.mat
% s contains data about singularities and the first and last points of a
% previous continuation run
path_to_this_script = get_path;
LP0_file = fullfile(path_to_this_script, 'Data', 'testbruss_LP0.mat');
load(LP0_file, 's');

%% Select the second singularity from the data loaded from testbruss_LP0.mat
% this singularity is a limit point
ID = 2;
data = s(ID).data;

%% Use L and A as the active parameters.
% The indices of L and A in the parameter array ( [N; L; A; B; Dx; Dy] ) are 2
% and 3. Hence the active parameter indices must be set to 2 and 3 to select L
% and A the active parameters.
active_parameter_indices = [2 3];


%% Initialize the continuation
% the init functions in cl_matcontL initialize the global variable cds, which is
% a struct contains information the continuer needs.
[x0,v0] = init_LP_LP_L(@brusselator_1d, [], [], active_parameter_indices, data);

%% Continuation
[singularities, datafile] = contL(@limitpointL,x0,v0,opt);

%% Plot results

x = loadPoint(datafile);

figure
hold on
title('testbruss\_LP1')
plot(x(end-1, :), x(end, :))

for i = 1:length(singularities)
  xii = singularities(i).data.x;
  plot(xii(end-1), xii(end), 'r.')
  text(xii(end-1), xii(end), singularities(i).label)
end
xlabel 'l'
ylabel 'a'   