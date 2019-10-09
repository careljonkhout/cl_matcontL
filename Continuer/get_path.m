function path = get_path
  stack         = dbstack('-completenames');
  path_elements = strsplit(stack(2).file, filesep);
  path          = strjoin(path_elements(1:end-1), filesep);
end