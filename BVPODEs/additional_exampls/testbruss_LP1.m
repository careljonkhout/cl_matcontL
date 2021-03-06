function testbruss_LP1()
% Test script for 1d brusselator (original MATCONT data)

% This example performs a continuation of limit points starting from a
% point located in testbruss_LP0

%% Options
opt = contset(); %Clear previous options
opt = contset(opt,'Cont_LogFile',                1);
opt = contset(opt,'Cont_DiagnosticsLevel',       3);
opt = contset(opt,'Cont_Direction',              0);
opt = contset(opt,'Cont_InitStepsize',           3);
opt = contset(opt,'Cont_MaxStepsize',            3);
opt = contset(opt,'Cont_MinStepsize',         1e-8);
opt = contset(opt,'Cont_MaxCorrIters',          12);
opt = contset(opt,'Cont_MaxNewtonIters',         5);
opt = contset(opt,'Cont_FunTolerance',        1e-5);
opt = contset(opt,'Cont_VarTolerance',        1e-5);
opt = contset(opt,'Cont_SmoothingAngle',     pi/30);
opt = contset(opt,'Cont_Singularities',          1);
opt = contset(opt,'Cont_Userfunctions',          0);
opt = contset(opt,'Cont_MaxNumPoints',         150);
opt = contset(opt,'Loc_UseLocators',       [1 1 1]);
opt = contset(opt,'Loc_Testf_MaxIters',         40); 
opt = contset(opt,'Loc_Testf_FunTolerance',   1e-4);
opt = contset(opt,'Loc_Testf_VarTolerance',   1e-4);
opt = contset(opt,'Loc_Testf_MaxIters',         16);
opt = contset(opt,'CIS_SparseSolvers',           1);
opt = contset(opt,'CIS_NStableRef',              9);
opt = contset(opt,'CIS_NExtra',                  4);
opt = contset(opt,'RIC_Cayley_Shift',           10);

opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename',     'testbruss_LP1');

%% Continuation
load('Data\testbruss_LP0')
ap = [2, 3];
ID = 2;
if(strcmp(s(ID).label ,'LP'))
    data = s(ID).data;
    %
    [x0,v0]      = init_LP_LP_L(@bruss_1d, [], [], ap, data);
    contL(@limitpointL,x0,v0,opt);
    
    
    %% Plot results
    x = loadPoint('Data\testbruss_LP1.dat');
    load('Data\testbruss_LP1.mat')
    
    plot(x(end-1, :), x(end, :))
    hold on
    for ii = 1:length(s)
        xii = s(ii).data.x;
        plot(xii(end-1), xii(end), 'r.')
        text(xii(end-1), xii(end), s(ii).label)
    end
end    

xlabel 'l'
ylabel 'a'