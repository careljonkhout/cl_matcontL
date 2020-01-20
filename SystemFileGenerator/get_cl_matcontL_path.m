function path = get_cl_matcontL_path()
  stack         = dbstack('-completenames');
  path_elements = strsplit(stack(1).file, filesep);
  path          = strjoin(path_elements(1:end-2), filesep);
end