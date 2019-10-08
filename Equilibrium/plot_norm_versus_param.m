% meant to be passed as a function handle to contL

function plot_norm_versus_param(currpoint, trialpoint)
  if ~ usejava('jvm')
    % When debugging mex files with gdb (GNU debugger) on linux, we ussually
    % start matlab without jvm (java virtual machine). This clause prevents a
    % continuation from raising an error in this method, when running without a
    % jvm (i.e. when debugging mex files with gdb on linux)
    return
  end
  if isempty(findobj('type', 'figure'))
    figure
  end
  hold on;
  curr_norm      = norm(currpoint.x(1:end-1));
  trial_norm     = norm(trialpoint.x(1:end-1));
  curr_param     = currpoint.x(end);
  trial_param    = trialpoint.x(end);
  

  plot([curr_param trial_param],[curr_norm trial_norm],'b');
  drawnow limitrate
end