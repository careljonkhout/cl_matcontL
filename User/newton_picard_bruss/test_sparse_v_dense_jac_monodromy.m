function test_sparse_v_dense_jac_monodromy
  clear global cds
  clear global lds
  N = 25;                     
  odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
  a = -1;
  b = -0.3;
  q_inf = -0.72;
  parameters = {a ; b; q_inf};
  load('/home/carel/Documents/cl_matcontL/User/fusion_newton_picard/Data/fusion_np_from_previous_20-Mar-2019_18_43_33/point 3.mat', 'point')
  x0 = point.x(1:end-2);
  period = point.x(end-1);
  parameters{3} = point.x(end);
  global contopts cds
  tol = 1e-11;
  contopts.abs_tol = tol;
  contopts.rel_tol = tol;
  contopts.monodromy_map_abs_tol = tol;
  contopts.monodromy_map_rel_tol = tol;
  
  cds.nphases = 3*(N-1);

  handles = feval(odefile);

  f =@(t, y) feval(handles{2},t, y, parameters{:});
  
  cds.integrator = @ode15s;
  cds.jacobian_ode = handles{3};
  
  pattern = feval(handles{3},0,ones(cds.nphases,1),parameters{:});
  
  integration_opt = odeset(...
    'AbsTol',      contopts.abs_tol,    ...
    'RelTol',      contopts.rel_tol,    ...
    'Jacobian',     @(t,y) feval(handles{3},t,y,parameters{:}) ...
  );
  
  cds.cycle_orbit  = ode15s(f, [0 period], x0, integration_opt);
  
  if false
  monodromy_map = @(x) monodromy_map_finite_differences( ...
    x0, x, period, parameters);
  disp('by finite differences')
  tic
  [~, multiplier_matrix] = eigs(monodromy_map, cds.nphases,5);
  toc
  multipliers = diag(multiplier_matrix);
  disp(multipliers);
  
  monodromy_map = @(x) NewtonPicard.SingleShooting.monodromy_map( ...
    x, period, parameters);
  
  disp('by integrating dense jacobian');
  tic
  monodromy_map(f(0,x0));
  toc


  end
  odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d_sparse_jac', N));
  handles = feval(odefile);
  
  i=1;
  monodromy_map = @(x) monodromy_map_from_jac(x, period, parameters);
  
  
  disp('by integrating sparse jacobian');
  tic
  Mx = monodromy_map(f(0,x0));
  toc
  plot(Mx ./ f(0,x0));

  pattern = feval(handles{3},0,ones(3*(N-1)),parameters{:})';
  monodromy_map = @(x) monodromy_map_pattern(x, period, parameters);
  disp('with JPattern');
  tic
  Mx = monodromy_map(f(0,x0));
  toc
  figure
  plot(Mx ./ f(0,x0));
  
  function x_end = shoot(x, period, ~)
    [~, orbit] = ode15s(f, [0 period], x, integration_opt);
    x_end = orbit(end,:)';
  end



  function d_phi__d_x = monodromy_map_finite_differences(...
      x_cycle, x, period, parameters)
    h = 1e-9;
    phi_1 = shoot(x_cycle+x*h, period, parameters);
    phi_2 = shoot(x_cycle-x*h, period, parameters);
    d_phi__d_x = (phi_2 - phi_1)/h/2; 

  end

  function Mx = monodromy_map_from_jac(phases_0, period, parameters)
    int_opt = odeset(...
      'AbsTol',       contopts.monodromy_map_abs_tol, ...
      'RelTol',       contopts.monodromy_map_rel_tol, ...
      'Jacobian',     @(t,y) feval(handles{3}, ...
                        t, deval(cds.cycle_orbit,t), parameters{:})' ...
    );                
    dydt_mon = @(t, y) ...
      feval(handles{3},t, deval(cds.cycle_orbit,t), parameters{:})' * y;
    tic
    [~,orbit] = ode15s(dydt_mon, [0 period], phases_0, int_opt);
    toc
    Mx = orbit(end,:)';
  end

  
  function Mx = monodromy_map_pattern(phases_0, period, parameters)
    int_opt = odeset(...
      'AbsTol',       contopts.monodromy_map_abs_tol, ...
      'RelTol',       contopts.monodromy_map_rel_tol, ...
      'JPattern',     pattern ...
    );                
    dydt_mon = @(t, y) ...
      feval(handles{3},t, deval(cds.cycle_orbit,t), parameters{:})' * y;
    [~,orbit] = ode15s(dydt_mon, [0 period], phases_0, int_opt);
    Mx = orbit(end,:)';
  end
end

