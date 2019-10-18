%% Continuation of a Hopf point
% This example performs a continuation of a Hopf point starting from a
% point located in testbruss_HP0.m (Chien 97 data)

function testbruss_HP1()

%% Options
opt = contset(    'contL_LogFile',           1);
opt = contset(opt,'contL_DiagnosticsLevel',  3);
opt = contset(opt,'Backward',               0);
opt = contset(opt,'InitStepsize',        1e-2);
opt = contset(opt,'MaxStepsize',         1e-2);
opt = contset(opt,'MinStepsize',         5e-8);
opt = contset(opt,'MaxCorrIters',          10); % default 10
opt = contset(opt,'MaxNewtonIters',         5);
opt = contset(opt,'FunTolerance',        1e-5);
opt = contset(opt,'VarTolerance',        1e-5);
opt = contset(opt,'contL_SmoothingAngle',pi/30);
opt = contset(opt,'Singularities',           1);
opt = contset(opt,'Userfunctions',           0);
opt = contset(opt,'MaxNumPoints',           70);
opt = contset(opt,'Locators',        [1 1 1 1]);
opt = contset(opt,'MaxTestIters',           10);
opt = contset(opt,'contL_Testf_FunTolerance', 1e-5);
opt = contset(opt,'contL_Testf_VarTolerance', 1e-4);
opt = contset(opt,'CIS_SparseSolvers',      1);
opt = contset(opt,'CIS_NStableRef',        6);
opt = contset(opt,'CIS_NExtra',             4);
opt = contset(opt,'CIS_Ric_Cayley_Shift',       10);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename', 'testbruss_HP1');

%% Load the struct s from the file Data/testbruss_HP0.mat
% s contains data about singularities and the first and last points of a
% previous continuation run
path_to_this_script = get_path;
HP0_file = fullfile(path_to_this_script, 'Data', 'testbruss_HP0.mat');
load(HP0_file, 's');

%% Select the second singularity from the data loaded from testbruss_U0.mat
% this singularity is a Hopf point
ID = 2;
data = s(ID).data;

%% Use Dx and Dy as the active parameters.
% The indices of Dx and Dy in the parameter array ( [N; L; A; B; Dx; Dy] ) are
% 4 and 5. Hence the active parameter indices must be set to 4 and 5 to select
% Dx and Dy as the active parameters.
active_parameter_indices = [4, 5];

%% Initialize the continuation
% the init functions in cl_matcontL initialize the global variable cds, which is
% a struct contains information the continuer needs.
[x0,v0] = init_H_H_L(@brusselator_1d, [], [], active_parameter_indices, data);

%% Continuation
[singularities, datafile] = contL(@hopfL,x0,v0,opt);


x = loadPoint(datafile);

figure
hold on
title('testbruss\_HP1')
plot(x(end-2, :), x(end-1, :));

for sing = singularities
  plot(x(end-2, sing.index), x(end-1, sing.index), 'r.');
  text(x(end-2, sing.index), x(end-1, sing.index), sing.label);
end

xlabel 'b'
ylabel 'Dx'