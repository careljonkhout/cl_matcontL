function testThyroid_7d1()

testThyroid_7d0
% This example performs continuation of equilibrium curve and locates 
% a limit point.

%% Options
opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',               1);
opt = contset(opt,'contL_DiagnosticsLevel',      3);
opt = contset(opt,'Backward',                    1);  % {0,1} = {Backward, Forward}
opt = contset(opt,'InitStepsize',             1e-2);   
opt = contset(opt,'MaxStepsize',              1e-1);  
opt = contset(opt,'MinStepsize',              1e-5);

opt = contset(opt,'MaxCorrIters',                5);  
opt = contset(opt,'MaxNewtonIters',              5);  
%opt = contset(opt,'FunTolerance',            1e-6); 
%opt = contset(opt,'VarTolerance',            1e-5); 
opt = contset(opt,'FunTolerance',             1e-8); 
opt = contset(opt,'VarTolerance',             1e-7);
opt = contset(opt,'MaxNumPoints',               10); 

opt = contset(opt,'Singularities',               1); 

opt = contset(opt,'Locators',              [1 0 0 0]);
opt = contset(opt,'MaxTestIters',                30); 
opt = contset(opt,'contL_Testf_FunTolerance',   1e-7); 
opt = contset(opt,'contL_Testf_VarTolerance',   1e-6); 
opt.Singularities
opt = contset(opt,'Multipliers',                 1);
opt = contset(opt,'Adapt',                       1);

%opt = contset(opt,'contL_EQ_BranchingMethod',       2);
%opt = contset(opt,'Userfunctions',               0); 
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename',     'testThyroid_7d1');

%% Continuation 


load('Data\testThyroid_7d0', 's')

ap = 12;
ID = 2;

if(strcmp(s(ID).label ,'H '))
    data = s(ID).data;
    x  = data.x;    
    % get the coordinate values of the hopf point:
    % from the continuation performed in testThyroid_7d0.m
    x1 = x(1:end-1);
    % get parameter values from start of previous continuation:
    p = data.P0; 
    % set active parameter to the value where the Hopf bifurcation
    % occurs:
    p(ap) = x(end);
    % limit cycle detection tolerance:
    h = 1e-6;
    % number of mesh points:
    ntst = 20;
    % number of colocation points per mesh interval:
    ncol = 4;
    [x0,v0] = init_H_LC_L(@Thyroid_7d0, x1, p, ap, h, ntst, ncol);
    [~, datafile] = contL(@limitcycleL,x0,v0,opt);
end

%% Plot results
figure
axes
[xlc, vlc, ~] = loadPoint(datafile); % DV: load computed cycles
load('Data\testThyroid_7d1', 's')            % DV: load singular points
plotcycle(xlc,vlc,s,[size(xlc,1) 1 2]);
pause

x = loadPoint('Data\testThyroid_7d0.dat');
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
