% recompiles all mex files

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

for i=1:length(source_files)
  compile(path, source_files{i})
end

function compile(path, file)
  compile_options = ['-outdir ' path];
  if ~isempty(regexp(mexext,'64','match'))
    compile_options = [compile_options ' -largeArrayDims'];
  end
  file_w_path = fullfile(path,file);
  compile_eval_string = ['mex ' compile_options ' -O ' file_w_path '.c '];
  fprintf(...
    ['The following command is being executed:\n' compile_eval_string '\n'])
  eval(compile_eval_string);
end
