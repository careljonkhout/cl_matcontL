%RUN ME FIRST!
if exist('contL', 'file')
  % in this case some other version of matcont might be on the path
  restoredefaultpath
  clearvars
end


%addpath(genpath(pwd))

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
addpath([cd '/LimitCycle/pqzschur']);
addpath([cd '/LimitCycleCodim2/']);
addpath([cd '/MultilinearForms/']);
addpath([cd '/Options/']);
addpath([cd '/PeriodDoubling/']);
addpath([cd '/SystemFileGenerator/']);
addpath([cd '/Systems/']);
addpath([cd '/TimeIntegration/']);



source_files = {
  'BVP_LC_jac';
  
  'BVP_PD_jac';
  
  'BVP_BPC_jacC';
  'BVP_BPC_jacCC';
  
  'BVP_LPC_jac';
  'BVP_LCX_jac';
  
  'BVP_NS_jac';
};

path = mfilename('fullpath');
path = path(1:end-length(mfilename));
path = fullfile(path,'LimitCycle');

for i=1:length(source_files)
  compile_if_needed(fullfile(path,source_files{i}))
end

function compile_if_needed(file)
  if ~(exist(strcat(file, '.', mexext),'file'))
    compile(file);
  end  
end

function compile(file)
  path = mfilename('fullpath');
  path = path(1:end-length(mfilename));
  path = fullfile(path,'LimitCycle');
  compile_options = ['-outdir ' path];
  if ~isempty(regexp(mexext,'64','match'))
    compile_options = [compile_options ' -largeArrayDims'];
  end
  compile_eval_string = ['mex ' compile_options ' -O ' file '.c '];
  fprintf(...
    ['The following command is being executed:\n' compile_eval_string '\n'])
  eval(compile_eval_string);
end
