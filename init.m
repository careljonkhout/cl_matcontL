%RUN ME FIRST!
function init

  %addpath(genpath(pwd))

  addpath([pwd '/BogdanovTakens/']);
  addpath([pwd '/BranchPointCycle/']);
  addpath([pwd '/CIS/']);
  addpath([pwd '/Continuer/']);
  addpath([pwd '/DataStorage/']);
  addpath([pwd '/Equilibrium/']);
  addpath([pwd '/Hopf/']);
  addpath([pwd '/LimitPoint/']);
  addpath([pwd '/Assertions/']);
  addpath([pwd '/LimitPointCycle/']);
  addpath([pwd '/LimitCycle/']);
  addpath([pwd '/LimitCycle/pqzschur']);
  addpath([pwd '/LimitCycleCodim2/']);
  addpath([pwd '/MultilinearForms/']);
  addpath([pwd '/Options/']);
  addpath([pwd '/PeriodDoubling/']);
  addpath([pwd '/SystemFileGenerator/']);
  addpath([pwd '/Systems/']);
  addpath([pwd '/primme/Matlab']);

  source_files = {
    'BVP_LC_jac';
    
    'BVP_PD_jac';
    
    'BVP_BPC_jacC';
    'BVP_BPC_jacCC';
    
    'BVP_LPC_jac';
    'BVP_LCX_jac';
    
    'BVP_NS_jac';
  };

  my_path = mfilename('fullpath');
  my_path = my_path(1:end-length(mfilename));
  my_path = fullfile(my_path,'LimitCycle');
  
  cd(my_path)

  for i=1:length(source_files)
    compile_if_needed(source_files{i})
  end
  
  my_path = mfilename('fullpath');
  my_path = my_path(1:end-length(mfilename));
  my_path = fullfile(my_path,'SystemFileGenerator');
  
  cd(my_path)
  
  compile_if_needed('replace_symbols');
  
  function compile_if_needed(file)
    if ~(exist(strcat(file, '.', mexext),'file'))
      compile(file);
    end  
  end

  function compile(file)
    % we quote the path to support spaces in the path
    compile_options = [];
    if ~isempty(regexp(mexext,'64','match'))
      compile_options = [compile_options ' -largeArrayDims'];
    end
    % we quote the path to support spaces in the path
    compile_eval_string = ['mex ' compile_options ' -O ''' file '.c'''];
    disp('The following command is being executed:');
    disp(compile_eval_string);
    eval(compile_eval_string);
  end

  cd ..
end
  
