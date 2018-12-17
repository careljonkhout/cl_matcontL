function testThyroid_7d0()
% Test script for 1d brusselator (original MATCONT data)

% This example performs continuation of equilibrium curve and locates 
% a limit point.

%% Options
opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',               1);
opt = contset(opt,'contL_DiagnosticsLevel',      3);
opt = contset(opt,'Backward',                    1);  % {0,1} = {Backward, Forward}
opt = contset(opt,'InitStepsize',             1e-5);   
opt = contset(opt,'MaxStepsize',              1e-5);  
opt = contset(opt,'MinStepsize',              1e-8);

opt = contset(opt,'MaxCorrIters',                5);  
opt = contset(opt,'MaxNewtonIters',              5);  
%opt = contset(opt,'FunTolerance',            1e-6); 
%opt = contset(opt,'VarTolerance',            1e-5); 
opt = contset(opt,'FunTolerance',             1e-8); 
opt = contset(opt,'VarTolerance',             1e-7);
opt = contset(opt,'MaxNumPoints',               102); 

opt = contset(opt,'CIS_SparseSolvers',          0);
%opt = contset(opt,'RIC_Cayley_Shift',              10);
opt = contset(opt,'CIS_NUnstable',             -1);
opt = contset(opt,'CIS_NExtra',                 0);
opt = contset(opt,'CIS_NSub',                   7); % NUnstable + NStableRef
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
opt = contset(opt, 'Filename',     'testThyroid_7d0');

%% Continuation 
N = 7;   
v0=1;
v1=1;
v01=1;
v2=1;  % 0 to 2 
v3=1;
v4=1;  %
v5=1;  % 
v6=1;  %
v7= 0;  %
b1 = 1;
a1 = 0.001;
Ss=100;
%Ni = 6;   % u(Ni) is a component of u to be monitored for homeostasis

%p = [N;v0;v1;v01;v2;v3;Ss];  ap1 = 2;
p = [N, v0, v1, v01, v2, v3, v4, v5, v6, v7, b1, a1, Ss];  
ap1 = 12;

[x0,v0]      = init_EP_EP_L(@Thyroid_7d0, [], p, ap1);
%contL(@equilibriumL_hom,x0,v0,opt);
[s,datafile] = contL(@equilibriumL,x0,v0,opt);

%% Plot results
x = loadPoint(datafile);
load(fullfile('Data','testThyroid_7d0.mat'),'s')
%N = s(1).data.P0(1);
%xx = s(2).data.x(:, 1)
%pause
plot(x(end, :), x(2, :)); % y2
hold on
plot(x(end, :), x(4, :)); % y4
for sii = s
    plot(x(end, sii.index), x(2, sii.index), 'r.');  % y2
    text(x(end, sii.index), x(2, sii.index), sii.label);  % y2
    plot(x(end, sii.index), x(4, sii.index), 'r.');  % y2
    text(x(end, sii.index), x(4, sii.index), sii.label);  % y2
end

xlabel 'v4'
ylabel 'y(2) = FT3'
%ylabel 'y(4) = TSH'
