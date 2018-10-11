function test_arch_BP1()
% Modeled on testbruss_BP1

opt = contset(); %Clear previous options
opt = contset(opt,'Cont_LogFile',          1);
opt = contset(opt,'Cont_DiagnosticsLevel', 3);
opt = contset(opt,'Backward',              1);
%%opt = contset(opt,'InitStepsize',  2e-1);
opt = contset(opt,'InitStepsize',  2e-1);  %d=3
%%opt = contset(opt,'InitStepsize',  1e-1);  %d=6
%%opt = contset(opt,'InitStepsize',  1e-1);  %d=12
%%opt = contset(opt,'MaxStepsize',   1e-0);
opt = contset(opt,'MaxStepsize',   3e-0);  %d=3
%%opt = contset(opt,'MaxStepsize',   6e-0);  %d=6, 12
opt = contset(opt,'MinStepsize',   1e-5);
opt = contset(opt,'MaxCorrIters',    12);
opt = contset(opt,'MaxNewtonIters',   6);
opt = contset(opt,'FunTolerance',  1e-6);
opt = contset(opt,'VarTolerance',  1e-5);
opt = contset(opt,'Cont_SmoothingAngle', pi/30);
opt = contset(opt,'Singularities',    1);
opt = contset(opt,'Userfunctions',    0);
%%opt = contset(opt,'MaxNumPoints',   500);
opt = contset(opt,'MaxNumPoints',   20);
opt = contset(opt,'Locator_MaxIters',          15);
opt = contset(opt,'Locator_FunTolerance',   1e-6);
opt = contset(opt,'Locator_VarTolerance',   1e-5);
%%opt = contset(opt,'Loc_VarTolerance',   1e-3);  %d=12
opt = contset(opt,'CIS_SparseSolvers',     1);
opt = contset(opt,'CIS_NStableRef',        6);
opt = contset(opt,'RIC_Cayley_Shift',     10);
opt = contset(opt,'EQ_BranchingMethod',    2);
opt = contset(opt,'Locators', [1 1 1]);
opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt, 'Filename', 'test_arch_BP1');

load('Data\test_arch_BP0')
ID = 2;
if strcmp(s(ID).label,'BP')
    data     = s(ID).data;
    [x0, v0] = init_BP_EP_L(@archcurve, [], [], [], data);
    contL(@equilibriumL, x0, v0, opt);
    
    
    global arch_mesh_setup;
    x = loadPoint('Data\test_arch_BP1.dat');
    load('Data\test_arch_BP1.mat');
    hold on
    plot(x(arch_mesh_setup.vdof, :), x(length(arch_mesh_setup.Iactive)+1, :));
    for sii = s
        plot(x(arch_mesh_setup.vdof, sii.index), x(length(arch_mesh_setup.Iactive)+1, sii.index), 'r.');
        text(x(arch_mesh_setup.vdof, sii.index), x(length(arch_mesh_setup.Iactive)+1, sii.index), sii.label);
    end
end

