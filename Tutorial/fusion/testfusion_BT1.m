function testfusion_BT1()
% Test script for the fusion model (original MATCONT data)

% This example performs continuation of the limit point located in
% testfusion_BT0 and locates a Bogdanov-Takens point.

%% Options
opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',          1);
opt = contset(opt,'contL_DiagnosticsLevel', 4);
opt = contset(opt,'Backward',              1);
opt = contset(opt,'InitStepsize',        0.1);
opt = contset(opt,'MaxStepsize',         0.5);
opt = contset(opt,'Increment',          1e-5);
opt = contset(opt,'FunTolerance',       1e-7);
opt = contset(opt,'VarTolerance',       1e-5);
opt = contset(opt,'contL_SmoothingAngle', pi/10);
opt = contset(opt,'Singularities',         1);
opt = contset(opt,'Userfunctions',         0);
opt = contset(opt,'MaxNumPoints',         10);
opt = contset(opt,'contL_ParallelComputing', 0);
opt = contset(opt,'CIS_SparseSolvers',     1);
opt = contset(opt,'CIS_NStableRef',        6); 
opt = contset(opt,'MaxTestIters',          40);  % more iterations for location
opt = contset(opt,'contL_Testf_FunTolerance',   1e-5);
opt = contset(opt, 'Filename', 'testfusion_BT1');  

%% load startpoint

path_to_this_script = get_path;
BT0_file = [path_to_this_script, '/Data/testfusion_BT0.mat'];
load(BT0_file, 's');
ap = [3, 8];
data = s(3).data;

%% Continuation

[x0,v0]                   = init_LP_LP_L(@fusion, [], [], ap, data);
[singularities, datafile] = contL(@limitpointL,x0,v0,opt);

%% Plot results

x = loadPoint(datafile);

figure
hold on
title('testfusion\_BT1')
plot(x(end, :), -x(end-1, :))

for i = 1:length(singularities)
  x_singularity = singularities(i).data.x;
  plot(x_singularity(end), -x_singularity(end-1), 'r.')
  text(x_singularity(end), -x_singularity(end-1), singularities(i).label)
end

xlabel 'b'
ylabel '-q_\infty'