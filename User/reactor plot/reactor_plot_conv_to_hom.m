

fig = figure;
hold on




for i=200:5:495
  
  file = sprintf('point %d.mat', i);
  path = fullfile(get_path(), '..', 'na_tubular_reactor', 'Data',  ...
    'na_tubular_reactor_orb_lc_15-Apr-2019_11_23_25', file);
  load(path)
  ncol = point.ncol;
  ntst = point.ntst;
  fine_mesh = get_fine_mesh(point.timemesh, ntst, ncol);
  pause(0.1)
  plot(fine_mesh, point.x(80 : 100 : end-2),'b');
end

xlabel('x_{80}')
ylabel('phi_{80}')

