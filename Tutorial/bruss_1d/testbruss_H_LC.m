function testbruss_HP0()
% Test script for 1d brusselator (Chien 97 data (!))

% This test function starts from a constant solution of the brussalator and
% detects a hopf point (and possibly branch points that are believed to be
% spurious solutions due to discretizaion.)

%% Options
opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',          1);
opt = contset(opt,'contL_DiagnosticsLevel', 3);
opt = contset(opt,'Backward',               0); %new
opt = contset(opt,'InitStepsize',         0.1);
opt = contset(opt,'MaxStepsize',          0.1);
opt = contset(opt,'MinStepsize',         1e-5);
opt = contset(opt,'Singularities',          1);
opt = contset(opt,'Userfunctions',          0);
opt = contset(opt,'MaxNumPoints',          10);
opt = contset(opt,'Locators',         [1 1 1]);
opt = contset(opt,'CIS_SparseSolvers',      1);
opt = contset(opt,'CIS_MaxUnstable',        5); %new
opt = contset(opt,'CIS_NStableRef',         6);
opt = contset(opt,'CIS_NExtra',             6);
opt = contset(opt,'CIS_Ric_Cayley_Shift',  10);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename', 'testbruss_HP0');

%% Start point
N = 15; L = 12; A = 4; B = 17.1; Dx = 1; Dy = 2;
p = [N; L; A; B; Dx; Dy];
ap = 4;

% construct starting point
x0 = zeros(2*N,1);
for i=1:N
  x0(i)   = A; 
  x0(N+i) = B/A;
end

%% Continuation
[x0,v0]      = init_EP_EP_L(@brusselator_1d, x0, p, ap);
contL(@equilibriumL,x0,v0,opt);

%% Continue limit cycle from Hopf

opt = contset(opt,'Singularities',          0);

x = loadPoint('Data\testbruss_HP0.dat');
load('Data\testbruss_HP0.mat');

system_handle = @bruss_1d;

data = s(2).data;
x  = data.x;    
x1 = x(1:end-1);

value_of_active_parameter_at_hopf = x(end,s(2).index);
p(ap) = value_of_active_parameter_at_hopf;

h = 1e-6;
ntst = 20;
ncol = 4;
%
[x0,v0] = init_H_LC_L(system_handle, x1, p, ap, h, ntst, ncol);
%[~, datafile] = contL(@limitcycleL,x0,v0,opt);

