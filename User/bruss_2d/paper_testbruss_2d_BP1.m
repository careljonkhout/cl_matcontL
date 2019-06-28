function paper_testbruss_2d_BP1()
% Test script for 2d brussaltor

opt = contset(); %Clear previous options
opt = contset(opt,'Cont_LogFile',          1);
opt = contset(opt,'Cont_DiagnosticsLevel', 3);
%%opt = contset(opt,'Backward',        0); %
opt = contset(opt,'Backward',        1); %
opt = contset(opt,'InitStepsize',   0.5); % n = 200,     N=10
opt = contset(opt,'MaxStepsize',    0.5); % n = 200,     N=10
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
%%opt = contset(opt,'MaxNumPoints',    30);
opt = contset(opt,'MaxNumPoints',    30);

opt = contset(opt,'CIS_SparseSolvers',     1);
opt = contset(opt,'CIS_NStableRef',        4);
opt = contset(opt,'CIS_NExtra',            4);

opt = contset(opt,'RIC_Cayley_Shift',     10);
opt = contset(opt,'Locators', [1 1 1]);
opt = contset(opt,'TestPath',mfilename('fullpath'));

opt = contset(opt,'Userfunctions',    1);
UserInfo{1}.name   = 'userf1';
UserInfo{1}.state  = 1;
UserInfo{1}.label  = 'u1';

UserInfo{2}.name   = 'userf2';
UserInfo{2}.state  = 0;
UserInfo{2}.label  = 'u2';

UserInfo{3}.name   = 'userf3';
UserInfo{3}.state  = 0;
UserInfo{3}.label  = 'u3';

UserInfo{4}.name   = 'userf4';
UserInfo{4}.state  = 0;
UserInfo{4}.label  = 'u4';

opt = contset(opt,'UserFuncInfo',UserInfo);

opt = contset(opt, 'Filename', 'paper_testbruss_2d_BP1');

%N = 200, 400 1600, 6400, 12800, 25600
%%N = 12800;
%%N = 6400;
%%N = 1600;
%%N = 400;
%%N = 10;
%%N = 50; %%n = 5,000
%%N = 100; %%n = 20,000

load('Data\paper_testbruss_2d_BP0')
ID = 2;
if strcmp(s(ID).label,'BP')
    data     = s(ID).data;
    [x0, v0] = init_BP_EP_L(@bruss_2d2, [], [], [], data);
    contL(@equilibriumL, x0, v0, opt);
    
    % plot result
    hold on
    x = loadPoint('Data\paper_testbruss_2d_BP1.dat');
    load('Data\paper_testbruss_2d_BP1.mat');
    N = s(1).data.P0(1);
    plot(x(2*N^2+1, :), x(1, :));
    for sii = s
        plot(x(2*N^2+1, sii.index), x(1, sii.index), 'r.');
        text(x(2*N^2+1, sii.index), x(1, sii.index), sii.label);
    end
    
    xlabel 'b'
    ylabel 'u(x=0,y=0)'
end


