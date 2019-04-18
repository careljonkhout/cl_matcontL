

fig = figure;
hold on
odefile = @nonadiabatic_tubular_reactor;
handles = feval(odefile);


for i=500:500
  
  file = sprintf('point %d.mat', i);
  path = fullfile(get_path(), '..', 'na_tubular_reactor', 'Data',  ...
    'na_tubular_reactor_orb_lc_15-Apr-2019_11_23_25', file);
  load(path)
  ncol = point.ncol;
  ntst = point.ntst;
  x = linspace(0,1,50);
  plot(x, point.x(8001:2:8100),'b');
  plot(x, point.x(8002:2:8100),'r');
  xlabel('x')
  ylabel('y,\phi')
  legend(['y   ';'\phi'], 'Location','east')
  
end

saveas(fig,'/home/carel/Documents/Master Thesis/reactor EP hom.svg');



