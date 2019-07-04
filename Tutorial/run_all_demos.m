% runs all .m files in the subdirectories of the Tutorials folder

% the order in which the files are run is the order in which the dir matlab
% command lists files. There are no guarantees on the order in which the
% files are run. 

pause off

if false || true
  try
    run_demos(dir(get_path()));
  catch exception
    pause on
    rethrow(exception)
  end
else
  run_demos(dir(get_path()))
end

pause on

function run_demos(directories)
  for i = 1:length(directories)
    directory = directories(i);
    if ~ directory.isdir
      continue
    end
    if strcmp(directory.name, '.') || strcmp(directory.name, '..')
      continue
    end

    my_path = [directory.folder '/' directory.name];
    files   = dir(my_path);

    for j = 1:length(files)
      filename = files(j).name;
      if length(filename) > 2 && strcmpi(filename(end-1:end), '.m')
        disp(' ')
        disp(repmat('+',2,80))
        disp(' ')
        disp(['     Now running:         ' directory.name '/' filename])
        disp(' ')
        disp(repmat('+',2,80))
        disp(' ')
        run([my_path '/' filename])
        drawnow
        
      end
    end
  end
end