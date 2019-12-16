function test_lorenz

  % generates the odefile of the Lorenz system
  % see https://en.wikipedia.org/wiki/Lorenz_system

  close all

  name = 'lorenz';
  pars = 'sigma r b';
  time = 't';
  max_ord = 1;
  equations = {
  'x'' = sigma * (-x + y)'
  'y'' = r*x - y - x*z'
  'z'' = -b*z + x*y'
  };

  lorenz_system = SystemFileGenerator.new( ...
                          name,pars,time,max_ord,equations,'cvode');
  lorenz_system.generate_file()

  sigma = 10;
  r=80;
  b = 8/3;

  tol = 1e-10;
  



  [t,y]= feval(@lorenz.cvode, ...
    't_values',                linspace(0,10,100), ...
    'initial_point',           100*ones(3,1), ...
    'ode_parameters',          [sigma r b], ...
    'abs_tol',                 tol, ...
    'rel_tol',                 tol);

  if usejava('jvm')
    figure
    plot(t, y);
  end
  point_on_limitcycle = y(end,:);



  [t,y] = feval(@lorenz.cvode, ...
    't_values',                linspace(0,3,200), ...
    'initial_point',           point_on_limitcycle, ...
    'cycle_detection',         true, ...
    'lower_bound_period',      1, ...
    'point_on_limitcycle',     point_on_limitcycle, ...
    'ode_parameters',          [sigma r b], ...
    'abs_tol',                 tol, ...
    'rel_tol',                 tol);


  if usejava('jvm')
    figure
    plot(t, y);
  end
  

  period = t(end);



  [~, eigenvalues, ~]  = eigs(@monodromy_map, 3, 3);
  global contopts;
  contopts = contset();
  
  disp(multipliers2str(diag(eigenvalues)));

  function Mx = monodromy_map(x)

    [~,~,Mx] = feval(@lorenz.cvode, ...
      't_values',                [0 period], ...
      'initial_point',           point_on_limitcycle, ...
      'sensitivity_vector',      x, ...
      'ode_parameters',          [sigma r b], ...
      'abs_tol',                 tol, ...
      'rel_tol',                 tol);
  end
end