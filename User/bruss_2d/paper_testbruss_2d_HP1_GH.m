function paper_testbruss_2d_HP1_GH()
% Test script for 2d brussaltor

opt = contset();
opt = contset(opt,'Cont_LogFile',          1);
opt = contset(opt,'Cont_DiagnosticsLevel', 3);
%%opt = contset(opt,'Backward',         0);
opt = contset(opt,'Backward',         0 );
opt = contset(opt,'InitStepsize',     1e-3);
opt = contset(opt,'MaxStepsize',      2);
opt = contset(opt,'MinStepsize',   1e-4);
% opt = contset(opt,'Cont_SmoothingAngle', 2*pi);
opt = contset(opt,'MaxCorrIters',    12);
opt = contset(opt,'MaxNewtonIters',  10);
%%opt = contset(opt,'FunTolerance',  1e-4);
%%opt = contset(opt,'VarTolerance',  1e-4); 
opt = contset(opt,'FunTolerance',  1e-7);
opt = contset(opt,'VarTolerance',  1e-5); 
opt = contset(opt,'Singularities',    1);
opt = contset(opt,'Userfunctions',    0);
opt = contset(opt,'MaxNumPoints',    30);
opt = contset(opt,'Adapt',        3);

opt = contset(opt,'MaxTestIters',   25);
opt = contset(opt,'Loc_Testf_FunTolerance',  1e-7);
%%opt = contset(opt,'Loc_TestTolerance', 1e-2);
opt = contset(opt,'Loc_Testf_VarTolerance',  1e-4);
%%opt = contset(opt,'Loc_UseLocators', [1 1 1 1]);  %%

opt = contset(opt,'CIS_SparseSolvers',      1);
opt = contset(opt,'CIS_NStableRef',         10);
opt = contset(opt,'CIS_NExtra',             4);

opt = contset(opt,'RIC_Cayley_Shift',       10);
opt = contset(opt,'TestPath',mfilename('fullpath'));

opt = contset(opt, 'Filename', 'paper_testbruss_2d_HP1_GH');
    
load('Data\paper_testbruss_2D_HP0')
ID = 2;
if(strcmp(s(ID).label ,'H '))
    data = s(ID).data;
    ap = [4, 6];
    %% Continuation
    [x0,v0]      = init_H_H_L(@bruss_2d2, [], [], ap, data);
    contL(@hopfL,x0,v0,opt);
    
    %% Plot results
    x = loadPoint('Data\paper_testbruss_2d_HP1_GH.dat');
    load('Data\paper_testbruss_2d_HP1_GH.mat')
    
    plot(x(end-2, :), x(end-1, :));
    hold on
    for sii = s
        plot(x(end-2, sii.index), x(end-1, sii.index), 'r.');
        text(x(end-2, sii.index), x(end-1, sii.index), sii.label);
    end
    
    xlabel 'b'
    ylabel 'Dy'
end 