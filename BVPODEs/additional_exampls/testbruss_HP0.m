function testbruss_HP0()
% Test script for 1d brusselator (Chien 97 data (!))

% This test function starts from a constant solution of the brussalator and
% detects a hopf point (and possibly branch points that are believed to be
% spurious solutions due to discretizaion.)

%% Options
opt = contset(); %Clear previous options
opt = contset(opt,'Cont_LogFile',          1);
opt = contset(opt,'Cont_DiagnosticsLevel', 1);
opt = contset(opt,'Cont_Direction',        0); %new
opt = contset(opt,'Cont_InitStepsize',   0.1);
opt = contset(opt,'Cont_MaxStepsize',    0.1);
opt = contset(opt,'Cont_MinStepsize',   1e-5);
opt = contset(opt,'Cont_Singularities',    1);
opt = contset(opt,'Cont_Userfunctions',    0);
opt = contset(opt,'Cont_MaxNumPoints',    10);
opt = contset(opt,'Loc_UseLocators', [1 1 1]);
opt = contset(opt,'CIS_SparseSolvers',     1);
opt = contset(opt,'CIS_NStableRef',        6);
opt = contset(opt,'CIS_NExtra',            6);
opt = contset(opt,'RIC_Cayley_Shift',     10);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename', 'testbruss_HP0');

%% Start point
N = 200; L = 12; A = 4; B = 17.1; Dx = 1; Dy = 2;
p = [N; L; A; B; Dx; Dy];
ap1 = 4;

% construct starting point
x0 = zeros(2*N,1);
for i=1:N
  x0(i)   = A; 
  x0(N+i) = B/A;
end

%% Continuation
[x0,v0]      = init_EP_EP_L(@bruss_1d, x0, p, ap1);
contL(@equilibriumL,x0,v0,opt);

%% Plot results

x = loadPoint('Data\testbruss_HP0.dat');
load('Data\testbruss_HP0.mat');
N = s(1).data.P0(1);
plot(x(2*N+1, :), x(1, :));
hold on
for sii = s
    plot(x(2*N+1, sii.index), x(1, sii.index), 'r.');
    text(x(2*N+1, sii.index), x(1, sii.index), sii.label);
end

xlabel 'b'
ylabel 'u(x=0)'