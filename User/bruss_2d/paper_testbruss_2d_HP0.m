function paper_testbruss_2d_HP0()
% Test script for 2d brussaltor

opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile', 3);
opt = contset(opt,'contL_DiagnosticsLevel', 3);
opt = contset(opt,'Backward',        1); % 
%%opt = contset(opt,'Backward',        0); % 
opt = contset(opt,'InitStepsize',   0.1); % n = 200,     N=10
opt = contset(opt,'MaxStepsize',    0.1); % n = 200,     N=10 
opt = contset(opt,'MinStepsize',   1e-5); % n = 200,     N=10 
%%opt = contset(opt,'InitStepsize',   0.1); % n = 800;   N=20
%%opt = contset(opt,'MaxStepsize',    0.1); % n = 800;   N=20
%%opt = contset(opt,'InitStepsize',   0.2); % n = 3200,  N=40 
%%opt = contset(opt,'MaxStepsize',    0.2); % n = 3200;  N=40
%%opt = contset(opt,'InitStepsize',   0.2); % n = 12800  N=80
%%opt = contset(opt,'MaxStepsize',    0.2); % n = 12800  N=80
opt = contset(opt,'MaxCorrIters',    12);
opt = contset(opt,'MaxNewtonIters',   6);  %% not 15
opt = contset(opt,'FunTolerance',  1e-7);   %-6
opt = contset(opt,'VarTolerance',  1e-5);   %-6
opt = contset(opt,'contL_SmoothingAngle', pi/30);
opt = contset(opt,'Singularities',    1);
opt = contset(opt,'Userfunctions',    0);
opt = contset(opt,'MaxNumPoints',    10);

opt = contset(opt,'MaxTestIters',     10);
opt = contset(opt,'contL_Testf_FunTolerance',1e-7);
opt = contset(opt,'contL_Testf_VarTolerance',1e-4);

opt = contset(opt,'CIS_SparseSolvers',     1);
opt = contset(opt,'CIS_NStableRef',        4);
opt = contset(opt,'CIS_NExtra',            4);

opt = contset(opt,'CIS_Ric_Cayley_Shift',     10);
opt = contset(opt,'Locators', [1 1 1]);
opt = contset(opt,'TestPath',mfilename('fullpath'));

opt = contset(opt, 'Filename', 'paper_testbruss_2d_HP0');

%%N = 10;  %%n = 200
%%N = 50; %%n = 5,000
N = 31; %%n = 20,000

%%L = 1; A = 1; B = 23; Dx = 1; Dy = 0.02;
L = 1; A =2; B = 5.45; Dx = 0.002; Dy = 0.004;
p = [N; L; A; B; Dx; Dy];
ap1 = 4;

n = N*N;
x0 = zeros(2*n,1);

for i=1:n
  x0(i)   = A;
  x0(n+i) = B/A;
end


[x0,v0]      = init_EP_EP_L(@bruss_2d2, x0, p, ap1);
contL(@equilibriumL,x0,v0,opt);

%% Plot results

x = loadPoint('Data\paper_testbruss_2d_HP0.dat');
load('Data\paper_testbruss_2d_HP0.mat');
N = s(1).data.P0(1);
plot(x(2*N^2+1, :), x(1, :));
hold on
for sii = s
    plot(x(2*N^2+1, sii.index), x(1, sii.index), 'r.');
    text(x(2*N^2+1, sii.index), x(1, sii.index), sii.label);
end

xlabel 'b'
ylabel 'u(x=0,y=0)'
