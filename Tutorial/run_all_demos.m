function run_all_demos(disable_pause_reset)

  % runs all .m files in the subdirectories of the Tutorials folder

  % the order in which the files are run is the order in which the dir matlab
  % command lists files. There are no guarantees on the order in which the
  % files are run. 

  my_pause off

  if nargin == 1 && disable_pause_reset
    run_demos(dir(get_path()));
  else
    try
      run_demos(dir(get_path()));
    catch exception
      my_pause on
      rethrow(exception)
    end
  end

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
      files   = flip(dir(my_path));

      for j = 1:length(files)
        filename = files(j).name;
        if length(filename) > 2 && strcmpi(filename(end-1:end), '.m')
          disp(' ')
          disp(repmat('+', 2, 80))
          disp(' ')
          disp(['     Now running:         ' directory.name '/' filename])
          disp(' ')
          disp(repmat('+', 2, 80))
          disp(' ')
          run([my_path '/' filename])
          drawnow
          % sometimes the repeated plot commands in all the demos fail to
          % produce output or write to the wrong plot, therefore we pause to
          % reduce the chance of this happening.
          pause(0.3)
        end
      end
    end
  end
end