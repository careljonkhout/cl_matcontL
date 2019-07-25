function testcusp2d_2()
% Test script for 1d brusselator (original MATCONT data)

% This example performs continuation of equilibrium curve and locates 
% a limit point.

%% Options
opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',               1);
opt = contset(opt,'contL_DiagnosticsLevel',      3);
opt = contset(opt,'Backward',                    1);  % {0,1} = {Backward, Forward}
opt = contset(opt,'InitStepsize',             2e-1);   
opt = contset(opt,'MaxStepsize',              2e-1);  
opt = contset(opt,'MinStepsize',              1e-8);

opt = contset(opt,'MaxCorrIters',                5);  
opt = contset(opt,'MaxNewtonIters',              5);  
%opt = contset(opt,'FunTolerance',            1e-6); 
%opt = contset(opt,'VarTolerance',            1e-5); 
opt = contset(opt,'FunTolerance',             1e-8); 
opt = contset(opt,'VarTolerance',             1e-7);
opt = contset(opt,'MaxNumPoints',               57); 

opt = contset(opt,'CIS_SparseSolvers',          0);
%opt = contset(opt,'RIC_Cayley_Shift',              10);
opt = contset(opt,'CIS_NUnstable',             -1);
opt = contset(opt,'CIS_NExtra',                 0);
opt = contset(opt,'CIS_NSub',                   2); % NUnstable + NStableRef
opt = contset(opt,'CIS_resizeable_flag',        0); % resize subspace MP 2018
opt = contset(opt,'CIS_Ric_SubspaceSelect', 'eig'); % (default 'ric') MP 2018 
opt = contset(opt,'CIS_DetectOverlap',          0);

opt = contset(opt,'Singularities',               1); 
%opt = contset(opt,'Locators',              []);
opt = contset(opt,'Locators',              [1 1 1]);
%opt = contset(opt,'Locators',          [1 1 1 1 0]);
opt = contset(opt,'contL_Loc_MaxCorrIters',      6);
opt = contset(opt,'contL_Loc_FunTolerance',   1e-4);
opt = contset(opt,'contL_Loc_VarTolerance',   1e-4');  
opt = contset(opt,'contL_SmoothingAngle',     pi/30);

%opt = contset(opt,'Locators',       []);
%opt = contset(opt,'MaxTestIters',                30); 
%opt = contset(opt,'contL_Testf_FunTolerance',   1e-7); 
%opt = contset(opt,'contL_Testf_VarTolerance',   1e-6); 

opt = contset(opt,'contL_EQ_BranchingMethod',       2);
opt = contset(opt,'Userfunctions',               0); 
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename',     'testcusp2d_2');

%% Continuation 
path_to_this_script = get_path;
testcusp2d_1_file = ...
               get_latest_singularity_file(path_to_this_script, 'testcusp2d_1');
load(testcusp2d_1_file, 's');
ap = [2, 3];
ID = 2;
if ~ strcmp(s(ID).label ,'LP')
  error('No limit point found')
end

data = s(ID).data;
%

problem_file = @cusp2d;

[x0,v0]      = init_LP_LP_L(problem_file, [], [], ap, data);
[singularities, datafile] = contL(@limitpointL,x0,v0,opt);


%% Plot results

x = loadPoint(datafile);
%N = s(1).data.P0(1);

figure
hold on
title('cusp2d_2')
plot(x(end-1, :), x(end, :));

for i = 1:length(singularities)
  xii = singularities(i).data.x;
  plot(xii(end-1), xii(end), 'r.')
  text(xii(end-1), xii(end), singularities(i).label)
end
xlabel 'l'
ylabel 'a'   

%N = s(1).data.P0(1);
%xx = s(2).data.x(:, 1)
%pause

xlabel 'h'
ylabel 'y(1) = x'
