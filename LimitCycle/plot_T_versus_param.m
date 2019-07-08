% meant to be passed as a function handle to contL

function plot_T_versus_param(currpoint, trialpoint)
  if isempty(findobj('type', 'figure'))
    figure
  end
  hold on;
  curr_T      = currpoint.x(end-1);
  trial_T     = trialpoint.x(end-1);
  curr_param  = currpoint.x(end);
  trial_param = trialpoint.x(end);
  
%  print_diag(2,'\nplotting: [%.3e %.3e]  [%.3e %.3e] \n\n', ...
%                [curr_param trial_param],[curr_T trial_T])
  plot([curr_param trial_param],[curr_T trial_T],'b');
  drawnow
end