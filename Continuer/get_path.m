function path = get_path
  stack = dbstack('-completenames');
  caller = stack(2);
  path = caller.file;
  % remove filename and .m
  path = path(1:end-length(caller.name)-2);
end