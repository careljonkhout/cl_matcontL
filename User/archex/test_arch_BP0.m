function test_arch_BP0()

% Modeled on testbruss_BP0
init
global arch_mesh_setup
arch_mesh_setup = [];

opt = contset(); %Clear previous options
opt = contset(opt,'Cont_LogFile',          1);
opt = contset(opt,'Cont_DiagnosticsLevel', 3);
opt = contset(opt,'Backward',             0);
%%opt = contset(opt,'InitStepsize',     0.1);
opt = contset(opt,'InitStepsize',     2e-0);  %d=3,6
%%opt = contset(opt,'InitStepsize',     8e-0);  %d=12
opt = contset(opt,'MinStepsize',   1e-5);  
%%opt = contset(opt,'MaxStepsize',      0.1);
opt = contset(opt,'MaxStepsize',      3e-0);  %d=3
%%opt = contset(opt,'MaxStepsize',      6e-0);  %d=6
%%opt = contset(opt,'Cont_MaxStepsize',      12e-0);  %d=12
opt = contset(opt,'MaxCorrIters',    12);
opt = contset(opt,'MaxNewtonIters',   5);
opt = contset(opt,'FunTolerance',  1e-6);
opt = contset(opt,'VarTolerance',  1e-5);
opt = contset(opt,'Cont_SmoothingAngle', pi/30);
opt = contset(opt,'Singularities',    1);
opt = contset(opt,'Userfunctions',    0);
%%opt = contset(opt,'MaxNumPoints',  500);
opt = contset(opt,'MaxNumPoints',    20);  %d=3,6, 12
opt = contset(opt,'Locator_MaxIters',         10);
opt = contset(opt,'Locator_VarTolerance',   1e-3);
opt = contset(opt,'Locator_FunTolerance',   1e-5);
opt = contset(opt,'CIS_SparseSolvers',     1);
opt = contset(opt,'CIS_NStableRef',        6);
opt = contset(opt,'RIC_Cayley_Shift',     10);
opt = contset(opt,'EQ_BranchingMethod',    2);
opt = contset(opt,'Locators', [1 1 1]);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename', 'test_arch_BP0');

density = 3;
%%density = 6;
%%density = 12;
p  = [0; density];
ap = 1;
[x0, v0] = init_EP_EP_L(@archcurve, [], p, ap);
contL(@equilibriumL, x0, v0, opt);


%  length(arch_mesh_setup.Iactive) % d=3 - N=1353; d=6 - N=5013

x = loadPoint('Data\test_arch_BP0.dat');
load('Data\test_arch_BP0.mat');
plot(x(arch_mesh_setup.vdof, :), x(length(arch_mesh_setup.Iactive)+1, :));
hold on
for sii = s
    plot(x(arch_mesh_setup.vdof, sii.index), x(length(arch_mesh_setup.Iactive)+1, sii.index), 'r.');
    text(x(arch_mesh_setup.vdof, sii.index), x(length(arch_mesh_setup.Iactive)+1, sii.index), sii.label);
end
