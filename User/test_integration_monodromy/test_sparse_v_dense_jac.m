function test_sparse_v_dense_jac
  % continuation of cycles cycles in fusion system
  clc
  clear global cds
  clear global lds
  N = 3;                     
  odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
  a = -1;
  b = -0.3;
  q_inf = -0.72;
  load('/home/carel/Documents/cl_matcontL/User/fusion_newton_picard/Data/fusion_np_from_previous_20-Mar-2019_18_43_33/point 3.mat', 'point')
  if N==25
    x0 = point.x(1:end-2);
  else
    x0 = ones(3*(N-1),1);
  end
  global contopts
  contopts.abs_tol = 1e-9;
  contopts.rel_tol = 1e-9;
  
  parameters = {a ; b; q_inf};
  handles = feval(odefile);

  int_opt = odeset( ...
    'AbsTol',      contopts.abs_tol,    ...
    'RelTol',      contopts.rel_tol    ...
  );



  t_end = 30;
  dydt = handles{2};
  f =@(t, y) dydt(t, y, parameters{:});
  disp("with dense Jac")
  tic
  sol = ode15s(f, [0 t_end], x0, int_opt);
  toc
  disp(sol.stats)
 
  int_opt = odeset( ...
    'AbsTol',      contopts.abs_tol,    ...
    'RelTol',      contopts.rel_tol,    ...
    'JPattern',    feval(handles{3},0,x0,parameters{:}) ...
  );


  dydt = handles{2};
  f =@(t, y) dydt(t, y, parameters{:});
  disp("seulexMex with Jpattern")
  tic
  [t,x] = seulexMex(f, [0 t_end], x0, int_opt);
  toc
 
    int_opt = odeset( ...
    'AbsTol',      contopts.abs_tol,    ...
    'RelTol',      contopts.rel_tol,    ...
    'Jacobian',    @(t,x) feval(handles{3},0,x,parameters{:}) ...
  );

  int_opt.RecomputeJACFactor = 0.01;
  dydt = handles{2};
  f =@(t, y) dydt(t, y, parameters{:});
  disp("seulexMex with Jacobian")
  tic
  [t,x] = seulexMex(f, [0 t_end], x0, int_opt);
  toc
 
  
  int_opt = odeset( ...
    'AbsTol',      contopts.abs_tol,    ...
    'RelTol',      contopts.rel_tol,    ...
    'JPattern',    feval(handles{3},0,x0,parameters{:}));
 %   'Jacobian',    @(t,x) feval(handles{3},0,x,parameters{:}) ...
  %);


 
  disp("with Jpattern and Jacobian")
  tic
  sol = ode15s(f, [0 t_end], x0, int_opt);
  toc
  disp(sol.stats)
  
    int_opt = odeset( ...
    'AbsTol',      contopts.abs_tol,    ...
    'RelTol',      contopts.rel_tol,    ...
    'JPattern',    feval(handles{3},0,x0,parameters{:}) ...
  );

  disp("with Jpattern")
  tic
  sol = ode15s(f, [0 t_end], x0, int_opt);
  toc
  disp(sol.stats)
  
  
  int_opt = odeset( ...
    'AbsTol',      contopts.abs_tol,    ...
    'RelTol',      contopts.rel_tol    ...
  );


 
  disp("without jac")
  tic
  sol = ode15s(f, [0 t_end], x0, int_opt);
  toc
  disp(sol.stats)
  
   f =@(t, y) dydt(t, y, parameters{:});
   
  if (N==25 || N==3)
    odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d_sparse_jac', N));
    handles = feval(odefile);

    int_opt = odeset( ...
      'AbsTol',      contopts.abs_tol,    ...
      'RelTol',      contopts.rel_tol,    ...
      'Jacobian',    @(t,x) feval(handles{3},0,x,parameters{:})' ...
    );

    disp("with sparse jac")
    tic
    sol = ode15s(f, [0 t_end], x0, int_opt);
    toc
    disp(sol.stats)
  end
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  function x_end = shoot(x, period, parameters)
    f =@(t, y) cds.dydt_ode(t, y, parameters{:});
    integration_opt = odeset(...
      'AbsTol',      contopts.shoot_abs_tol,    ...
      'RelTol',      contopts.shoot_rel_tol,    ...
      'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
    );
    [~, orbit] = ode15s(f, [0 period], x, integration_opt);
    x_end = orbit(end,:)';
  end



  function d_phi__d_x = monodromy_map_finite_differences(...
      x_cycle, x, period, parameters)
    h = 1e-6;
    phi_1 = shoot(x_cycle+x*h, period, parameters);
    phi_2 = shoot(x_cycle-x*h, period, parameters);
    d_phi__d_x = (phi_2 - phi_1)/h/2; 

  end
end

