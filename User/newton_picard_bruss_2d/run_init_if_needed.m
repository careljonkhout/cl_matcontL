function run_init_if_needed
  if ~ exist('contL', 'file')
    fullpath = mfilename('fullpath');
    file_directory = fullpath(1:end-length(mfilename));
    cd([file_directory '../..']);
    init
    cd(file_directory)
  end