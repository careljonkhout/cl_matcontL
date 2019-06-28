function paper_testbruss_2d_BP5()
% Test script for 2d brussaltor

opt = contset(); %Clear previous options
opt = contset(opt,'Cont_LogFile',          1);
opt = contset(opt,'Cont_DiagnosticsLevel', 3);
%%opt = contset(opt,'Backward',          0); %
opt = contset(opt,'Backward',        1); %
opt = contset(opt,'InitStepsize',   2e-2); % n = 200,     N=10
opt = contset(opt,'MaxStepsize',    2e-2); % n = 200,     N=10 
opt = contset(opt,'MinStepsize',   1e-5); % n = 200,     N=10 
%%opt = contset(opt,'InitStepsize',   0.1); % n = 800;   N=20
%%opt = contset(opt,'MaxStepsize',    0.1); % n = 800;   N=20
%%opt = contset(opt,'InitStepsize',   0.2); % N = 3200,  N=40 
%%opt = contset(opt,'MaxStepsize',    0.2); % N = 3200;  N=40
%%opt = contset(opt,'InitStepsize',   0.2); % N = 12800  N=80
%%opt = contset(opt,'MaxStepsize',    0.2); % N = 12800  N=80
opt = contset(opt,'MaxCorrIters',    12);
opt = contset(opt,'MaxNewtonIters',   5);  %% not 15
opt = contset(opt,'FunTolerance',  1e-7);   %-6
%%opt = contset(opt,'VarTolerance',  1e-5);   %-6
opt = contset(opt,'VarTolerance',  1e-5);   %-6
opt = contset(opt,'Cont_SmoothingAngle', pi/30);
opt = contset(opt,'Singularities',    1);
opt = contset(opt,'Userfunctions',    0);
opt = contset(opt,'MaxNumPoints',    20);

opt = contset(opt,'CIS_SparseSolvers',     1);
opt = contset(opt,'CIS_NStableRef',        4);
opt = contset(opt,'CIS_NExtra',            4);

opt = contset(opt,'RIC_Cayley_Shift',     10);
opt = contset(opt,'Locators', [1 1 1]);
opt = contset(opt,'TestPath',mfilename('fullpath'));

opt = contset(opt, 'Filename', 'paper_testbruss_2d_BP5');

%N = 200, 400 1600, 6400, 12800, 25600
%%N = 12800;
%%N = 6400;
%%N = 1600;
%%N = 400;
%%N = 10;
%%N = 50; %%n = 5,000
%%N = 100; %%n = 20,000

load('Data\paper_testbruss_2d_BP3')
ID = 3; %user point
if strcmp(s(ID).label,'u3')
    data = s(ID).data;
    ap1  = 4;
    u = data.x(1:end-1);
    [x0, v0] = init_EP_EP_L(@bruss_2d2, u, data.P0, ap1);
    contL(@equilibriumL, x0, v0, opt);
    
    % plot result
    hold on
    x = loadPoint('Data\paper_testbruss_2d_BP2.dat');
    load('Data\paper_testbruss_2d_BP2.mat');
    N = s(1).data.P0(1);
    plot(x(2*N^2+1, :), x(1, :));
    for sii = s
        plot(x(2*N^2+1, sii.index), x(1, sii.index), 'r.');
        text(x(2*N^2+1, sii.index), x(1, sii.index), sii.label);
    end
    
    xlabel 'b'
    ylabel 'u(x=0,y=0)'
end