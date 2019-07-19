function path = get_path
  stack = dbstack('-completenames');
  caller = stack(2);
  fullpath = caller.file;
  
  stack = dbstack();
  caller = stack(2);
  filename = caller.file;
  % remove filename and .m
  path = fullpath(1:end-length(filename));
end