% Meant to be passed as a function handle to contL. This causes the amplitudes
% of all coordinates will be plotted, during single_shooting cycle
% continuations.

function plot_amplitudes_versus_param(currpoint, trailpoint)
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
  for i = 1 : cds.n_phases
    curr_amplitude      = currpoint.amplitudes(i);
    trial_amplitude     = trailpoint.amplitudes(i);
    curr_param          = currpoint.x(end);
    trial_param         = trailpoint.x(end);
  
    %print_diag(3,'\nplotting: [%.3e %.3e]  [%.3e %.3e] \n\n', ...
    %              [curr_param trial_param],[curr_T trial_T])
    plot([curr_param trial_param],[curr_amplitude trial_amplitude],'b');
  end
  drawnow limitrate
end