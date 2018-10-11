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
opt = contset(opt,'contL_ParallelComputing', 1);
opt = contset(opt,'CIS_SparseSolvers',     1);
opt = contset(opt,'CIS_NStableRef',        6); 
opt = contset(opt,'MaxTestIters',          40);  % more iterations for location
opt = contset(opt,'contL_Testf_FunTolerance',   1e-5);
opt = contset(opt, 'Filename', 'testfusion_BT1');  

%% load startpoint

load('Data\testfusion_BT0.mat') 
ap = [3, 8];
data = s(3).data;

%% Continuation

[x0,v0]      = init_LP_LP_L(@fusion, [], [], ap, data);
contL(@limitpointL,x0,v0,opt);

%% Plot results

x = loadPoint('Data\testfusion_BT1.dat');
load('Data\testfusion_BT1.mat') 
plot(x(end, :), -x(end-1, :))
hold on
for ii = 1:length(s)
    xii = s(ii).data.x;
    plot(xii(end), -xii(end-1), 'r.')
    text(xii(end), -xii(end-1), s(ii).label)
end

xlabel 'b'
ylabel '-q_\infty'