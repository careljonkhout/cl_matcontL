function testfusion_BT2()
% Test script for the fusion model (original MATCONT data)

% This example performs continuation of the Bogdanov-Takens point located
% in testfusion_BT1 and locates degeneracies along the curve

opt = contset(); %Clear previous options
opt = contset(opt,'contL_LogFile',          4);
opt = contset(opt,'contL_DiagnosticsLevel', 4);
opt = contset(opt,'Backward',              0);
opt = contset(opt,'InitStepsize',       0.05);
opt = contset(opt,'MaxStepsize',        0.05);
opt = contset(opt,'Increment',          1e-5);
opt = contset(opt,'FunTolerance',       1e-7);
opt = contset(opt,'VarTolerance',       1e-5);
opt = contset(opt,'contL_SmoothingAngle', pi/10);
opt = contset(opt,'Singularities',         1);
opt = contset(opt,'Userfunctions',         0);
opt = contset(opt,'MaxNumPoints',         30);
opt = contset(opt,'contL_ParallelComputing',0);
opt = contset(opt,'CIS_SparseSolvers',     1);
opt = contset(opt,'CIS_NStableRef',        6);
opt = contset(opt,'MaxTestIters',         30);  % more iterations for location
opt = contset(opt,'contL_Testf_FunTolerance',   1e-5);
opt = contset(opt, 'Filename', 'testfusion_BT2');  

%% load startpoint
path_to_this_script = get_path;
BT1_file = [path_to_this_script, 'Data/testfusion_BT1.mat'];
load(BT1_file, 's');
ap = [3, 7, 8];
data = s(2).data;

%% continuation

[x0,v0]                   = init_BT_BT_L(@fusion, [], [], ap, data);
[singularities, datafile] = contL(@bogdanovtakensL,x0,v0,opt);

%% plot


x = loadPoint(datafile);

plot3(x(end, :), -x(end-2, :), x(end-1, :))
hold on
for i = 1:length(singularities)
  x_sing = singularities(i).data.x;
  plot3(x_sing(end), -x_sing(end-2), x_sing(end-1), 'r.')
  text( x_sing(end), -x_sing(end-2), x_sing(end-1), singularities(i).label)
end

xlabel 'a'
ylabel '-q_\infty'
zlabel 'b'