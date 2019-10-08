function plot_singularity_of_equilibria(s)
  persistent vertical_alignment
  if ~ usejava('jvm')
    % When debugging mex files with gdb (GNU debugger) on linux, we ussually
    % start matlab without jvm (java virtual machine). This clause prevents a
    % continuation from raising an error in this method, when running without a
    % jvm (i.e. when debugging mex files with gdb on linux)
    return
  end
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
  norm_of_x         = norm(s.data.x(1:end-1));
  plot(parameter_value, norm_of_x, 'r*')
  text(parameter_value, norm_of_x, s.label, ...
        'VerticalAlignment', vertical_alignment)
	drawnow limitrate
end

