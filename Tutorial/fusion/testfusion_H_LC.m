
% Test script for the fusion model (original MATCONT data)

% This example performs continuation of equilibrium curve and locates a
% limit point.

%% Options


N         = 30;


opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',          1);
opt = contset(opt,'contL_DiagnosticsLevel', 4);
opt = contset(opt,'Backward',              0);
opt = contset(opt,'InitStepsize',       0.03);
opt = contset(opt,'MinStepsize',        1e-5);
opt = contset(opt,'MaxStepsize',        0.03);
opt = contset(opt,'Increment',          1e-5);
opt = contset(opt,'FunTolerance',       1e-5);
opt = contset(opt,'VarTolerance',       1e-4);
opt = contset(opt,'contL_SmoothingAngle',  pi/10);
opt = contset(opt,'Singularities',         1);
opt = contset(opt,'Userfunctions',         0);
opt = contset(opt,'MaxNumPoints',         20);
opt = contset(opt,'contL_ParallelComputing', 0);  % using parallel computing toolbox
opt = contset(opt,'CIS_SparseSolvers',     1);
opt = contset(opt,'CIS_NStableRef',        6);
opt = contset(opt,'CIS_Ric_Transform',   'cayley');
opt = contset(opt,'CIS_Ric_Cayley_Shift',      4);
opt = contset(opt,'contL_EQ_BranchingMethod',    2);
opt = contset(opt,'Locators',        [1 1 1]);
opt = contset(opt,'contL_Testf_FunTolerance', 1e-4);
opt = contset(opt, 'Filename', ['testfusion_BT0_sage');

%% Parameters

q_inf     = -0.008;

a         = -0.01;
b         = -0.01; 


p = [ a, b, q_inf];
ap1 = 3;  % active parameter

%% Continuation
filename = mfilename;
fullpath = mfilename('fullpath');
path_wo_filename = fullpath(1:end-length(filename));
fusion_systems_path = fullfile(path_wo_filename, ...
  '..', '..', 'Systems', 'fusion');
addpath(fusion_systems_path);

evalstr = ['system_handle = @fusion_precomputed_with_sage_N_' num2str(N)];
eval(evalstr);



%% Plot results
load(fullfile('Data',opt.Filename))

ap = 3;
data = s(2).data;
x  = data.x;    
x1 = x(1:end-1);
P0 = data.P0;
p  = [a;b;P0(2)];

h = 1e-6;
ntst = 20;
ncol = 4;
%
[x0,v0] = init_H_LC_L(system_handle, x1, p, ap, h, ntst, ncol);
[~, datafile] = contL(@limitcycleL,x0,v0,opt);