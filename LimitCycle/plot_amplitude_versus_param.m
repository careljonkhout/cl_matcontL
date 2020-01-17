% The result of plot_amplitude_versus_param is a function handle meant to be
% passed to contL. This causes the amplitude of the i-th coordinate to be
% plotted during single_shooting cycle continuations.


function func = plot_amplitude_versus_param(i)
  func = @(currpoint, trailpoint) ...
                plot_amplitude_i_versus_param(currpoint, trailpoint, i);
end

function plot_amplitude_i_versus_param(currpoint, trailpoint, i)
  if ~ usejava('jvm')
    % When debugging mex files with gdb (GNU debugger) on linux, we ussually
    % start matlab without jvm (java virtual machine). This clause prevents a
    % continuation from raising an error in this method, when running without a
    % jvm (i.e. when debugging mex files with gdb on linux)
    return
  end

  global cds

  switch func2str(cds.curve)
    case 'single_shooting'

    otherwise
      error('plot_max_amplitude_versus_param for %s is not yet implemented', ...
            func2str(cds.curve));
  end
  if isempty(findobj('type', 'figure'))
    figure
  end
  hold on;
  
  curr_amplitude      = currpoint.amplitudes(i);
  trial_amplitude     = trailpoint.amplitudes(i);
  curr_param          = currpoint.x(end);
  trial_param         = trailpoint.x(end);

  %print_diag(3,'\nplotting: [%.3e %.3e]  [%.3e %.3e] \n\n', ...
  %              [curr_param trial_param],[curr_T trial_T])
  plot([curr_param trial_param],[curr_amplitude trial_amplitude],'b');

  drawnow limitrate
end
