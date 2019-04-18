

fig = figure('Position',[0 0 800 600]);
hold on
odefile = @nonadiabatic_tubular_reactor;
handles = feval(odefile);


for i=648:648
  
  file = sprintf('point %d.mat', i);
  path = fullfile(get_path(), '..', 'na_tubular_reactor', 'Data',  ...
    'na_tubular_reactor_orb_lc_15-Apr-2019_11_23_25', file);
  load(path)
  ncol = point.ncol;
  ntst = point.ntst;
  fine_mesh = point.T * get_fine_mesh(point.timemesh, ntst, ncol);
  norm_of_gradient = zeros(size(fine_mesh));
  x_idx = 1:100;
  pause(0.2)
  p = num2cell(point.parametervalues);
  for j=1:length(norm_of_gradient)
    norm_of_gradient(j) = log10(norm(feval(handles{2},0,point.x(x_idx),p{:})));
    x_idx = x_idx + 100;  
  end
  plot(fine_mesh,norm_of_gradient,'b');
end
xlim([0 203])

xlabel('t')
ylabel('log_{10} of the norm of the gradient of the discretized system')

saveas(fig,'/home/carel/Documents/Master Thesis/reactor grad v t.svg');

