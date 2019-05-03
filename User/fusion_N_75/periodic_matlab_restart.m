nSessions = 100;
% test of limit point bifurcation detection using orhtogonal collocation

%fix computation time addition

% in the next line we assume that cl_matcontL init script init.m is located in
% the subdirectory of the directory in which this script is located.
cl_matcontL_path   = fullfile(get_path(),'..','..');
script_path        = get_path();
script             = 'extend_fusion_oc';


for i=1:nSessions
  exec_string = sprintf([
    % To avoid opening a window, and to minimize memory use we start matlab with
    % the option -nojvm. Note that in this mode, plotting is disabled, and any
    % command related to plotting (plot, xlim, figure, etc.) will immediatly
    % cause an error and stop the program.
    %
    % With the -nodisplay option, plot commands apparently do not cause errors,
    % but plots are not displayed.
    %
    % With the -nodesktop option, plots are shown.
    'matlab -nojvm -r "', ...
      'try;'...
        'cd /home/carel/Documents/cl_matcontL;', ...
        'init;', ...
        sprintf('cd ''%s'';', script_path), ...
        script ';', ...
      'catch e;', ...
        'disp(e.message);', ...
        'disp(e.stack);', ...
        'quit(1);', ...
      'end;', ...
      'quit;"']);
  errorlevel = system(exec_string);
  if errorlevel >= 1
    break
  end
end