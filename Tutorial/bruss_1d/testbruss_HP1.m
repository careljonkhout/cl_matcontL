function testbruss_HP1()
% Test script for 1d brusselator (Chien 97 data (!))

% This example performs a continuation of hopf points starting from a
% point located in testbruss_HP0.m

%% Options
opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',           1);
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

%% Find start point
ap = [4, 5];
load('Data/testbruss_HP0')

ID = 2;
if(strcmp(s(ID).label ,'H '))
    data = s(ID).data;
    
    %% Continuation
    [x0,v0]      = init_H_H_L(@brusselator_1d, [], [], ap, data);
    contL(@hopfL,x0,v0,opt);
    
    %% Plot results
    x = loadPoint('Data\testbruss_HP1.dat');
    load('Data\testbruss_HP1.mat')
    
    plot(x(end-2, :), x(end-1, :));
    hold on
    for sii = s
        plot(x(end-2, sii.index), x(end-1, sii.index), 'r.');
        text(x(end-2, sii.index), x(end-1, sii.index), sii.label);
    end
    
    xlabel 'b'
    ylabel 'Dx'
end