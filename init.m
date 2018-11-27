%RUN ME FIRST!
%restoredefaultpath
%clearvars

addpath([cd '/BogdanovTakens/']);
addpath([cd '/BranchPointCycle/']);
addpath([cd '/CIS/']);
addpath([cd '/Continuer/']);
addpath([cd '/DataStorage/']);
addpath([cd '/Equilibrium/']);
addpath([cd '/Hopf/']);
addpath([cd '/LimitPoint/']);
addpath([cd '/LimitPointCycle/']);
addpath([cd '/LimitCycle/']);
addpath([cd '/MultilinearForms/']);
addpath([cd '/Options/']);

addpath([cd '/SystemFileGenerator/']);
addpath([cd '/Systems/']);
addpath([cd '/TimeIntegration/']);

p = mfilename('fullpath');
p = p(1:length(p)-length(mfilename));
p = strcat(p,'/LimitCycle');
curdir = cd;
cd(p);

% Compile the c-files (optimized)
if ~(exist(strcat('BVP_LC_jac.',mexext),'file'))
if ~isempty(regexp(mexext,'64','match'))
    mex -largeArrayDims -O BVP_LC_jac.c;
    mex -largeArrayDims -O BVP_PD_jac.c;
    mex -largeArrayDims -O BVP_BPC_jacC.c;
    mex -largeArrayDims -O BVP_BPC_jacCC.c;
    mex -largeArrayDims -O BVP_LPC_jac.c;
    mex -largeArrayDims -O BVP_NS_jac.c;
    mex -largeArrayDims -O BVP_LCX_jac.c;
else
    mex -O BVP_LC_jac.c;
    mex -O BVP_PD_jac.c;
    mex -O BVP_BPC_jacC.c;
    mex -O BVP_BPC_jacCC.c;
    mex -O BVP_LPC_jac.c;
    mex -O BVP_NS_jac.c;
    mex -O BVP_LCX_jac.c;
end
end


% Return to directory we started in
cd (curdir);
