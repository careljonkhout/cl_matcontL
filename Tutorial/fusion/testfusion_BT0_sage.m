
% Test script for the fusion model (original MATCONT data)

% This example performs continuation of equilibrium curve and locates a
% limit point.

%% Options
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

%% Parameters

N=20;

q_inf     = -0.008;
a         = -0.01;
b         = -0.01; 

ap1 = 3;  % active parameter

datafile = [mfilename ...
  sprintf('_N_%d__a_%d__b_%d__q_inf_%d__ap_%d', N, a, b, q_inf, ap)];
datafile = strrep(datafile,'.','p');
% because matlab load function can't handle periods in filenames
disp('Output files are:');
disp([datafile '.dat']);
disp([datafile '.mat']);
opt = contset(opt, 'Filename', datafile);


p = [ a, b, q_inf];

%% Continuation
add_fusion_systems_to_path

evalstr = ['system_handle = @fusion_precomputed_with_sage_N_' num2str(N)];
eval(evalstr);

[x0,v0]      = init_EP_EP_L(system_handle, [], p, ap1); % initialization
contL(@equilibriumL,x0,v0,opt);                   % continuation

%% Plot results
x = loadPoint(fullfile('Data',[opt.Filename, '.dat']));
load(fullfile('Data',opt.Filename))
