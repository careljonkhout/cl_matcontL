function testcusp2d_1()
% Test script for 1d ... 

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
opt = contset(opt, 'Filename',     'testcusp2d_1');

%% Continuation 
N = 2;   
h = 0; r = 1;
p = [N, h, r];  
ap1 = 2;

problem_file = @cusp2d;
[x0,v0]      = init_EP_EP_L(problem_file, [], p, ap1);
[singularities, datafile] = contL(@equilibriumL,x0,v0,opt);

%% Plot results

x = loadPoint(datafile);
N = singularities(1).data.P0(1);

figure
hold on
title('cusp2d_1')
plot(x(N+1, :), x(1, :));

for singularity = singularities
    plot(x(N+1, singularity.index), x(1, singularity.index), 'r.');
    text(x(N+1, singularity.index), x(1, singularity.index), singularity.label);
end

xlabel 'h'
ylabel 'y(1) = x'
