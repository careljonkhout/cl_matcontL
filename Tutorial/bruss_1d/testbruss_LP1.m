%% Continuation of a limit point
% This example performs a continuation of a limit point starting from a point
% located in testbruss_LP0

function testbruss_LP1()


%% Options
opt = contset(    'contL_LogFile',               1);
%opt = contset(opt,'contL_DiagnosticsLevel',     3);
opt = contset(opt,'contL_DiagnosticsLevel',      5);
opt = contset(opt,'Backward',                    0);
opt = contset(opt,'InitStepsize',                3);
opt = contset(opt,'MaxStepsize',                 3);
opt = contset(opt,'MinStepsize',              1e-8);
opt = contset(opt,'MaxCorrIters',               12);
opt = contset(opt,'MaxNewtonIters',              5);
opt = contset(opt,'FunTolerance',             1e-5);
opt = contset(opt,'VarTolerance',             1e-5);
opt = contset(opt,'contL_SmoothingAngle',    pi/30);
opt = contset(opt,'Singularities',               1);
opt = contset(opt,'Userfunctions',               0);
opt = contset(opt,'MaxNumPoints',              150);
opt = contset(opt,'Locators',              [1 1 1]);
opt = contset(opt,'contL_Testf_FunTolerance', 1e-4);
opt = contset(opt,'contL_Testf_VarTolerance', 1e-4);
opt = contset(opt,'MaxTestIters',               16);
opt = contset(opt,'CIS_SparseSolvers',           1);
opt = contset(opt,'CIS_NStableRef',              9);
opt = contset(opt,'CIS_MaxUnstable',             5); % MP new
opt = contset(opt,'CIS_NExtra',                  4);
opt = contset(opt,'CIS_Ric_Cayley_Shift',       10);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename',     'testbruss_LP1');

%% Load the struct s from the file Data/testbruss_LP0.mat
% s contains data about singularities and the first and last points of a
% previous continuation run
path_to_this_script = get_path;
LP0_file = fullfile(path_to_this_script, 'Data', 'testbruss_LP0.mat');
load(LP0_file, 's');

%% Select the second singularity from the data loaded from testbruss_U0.mat
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