function plot_cycle_collocation(point, coords_to_plot)
  % number of "phase variables" (dependent variables of the system of ODEs):
  nphase = numel(point.multipliers); 
  ntst = point.ntst; % number of mesh intervals
  ncol = point.ncol; % number of collocation points
  fine_mesh = get_fine_mesh(point.timemesh, ntst, ncol);
  x = reshape(point.x(1:end-2), nphase, ntst * ncol + 1)';
  if nargin > 1
    x = x(:, coords_to_plot);
  end
  plot(fine_mesh * point.T, x - x(1,:));
end
  