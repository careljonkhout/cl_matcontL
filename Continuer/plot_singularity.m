function plot_singularity(s)
  persistent vertical_alignment
  if isempty(vertical_alignment)
    vertical_alignment = 'bottom';
    % switch vertical alignment of label after each bifurcation to minimize
    % overlap
  elseif strcmp(vertical_alignment, 'bottom')
    vertical_alignment = 'top';
  else
    vertical_alignment = 'bottom';
  end
  parameter_value = s.data.x(end);
  period          = s.data.x(end-1);
  plot(parameter_value, period, 'r*')
  text(parameter_value, period, s.label, ...
        'VerticalAlignment', vertical_alignment)
end

