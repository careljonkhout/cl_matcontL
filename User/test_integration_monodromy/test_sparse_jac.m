function test_sparse_v_dense_jac
  % continuation of cycles cycles in fusion system
  clc
  clear global cds
  clear global lds
  N = 25;                     
  odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
  a = -1;
  b = -0.3;
  q_inf = -0.72;
  load('/home/carel/Documents/cl_matcontL/User/fusion_newton_picard/Data/fusion_np_from_previous_20-Mar-2019_18_43_33/point 3.mat', 'point')
  x0 = point.x(1:end-2);
  global contopts
  contopts.abs_tol = 1e-13;
  contopts.rel_tol = 1e-13;
  
  parameters = {a ; b; q_inf};
  handles = feval(odefile);

  int_opt = odeset( ...
    'AbsTol',      contopts.abs_tol,    ...
    'RelTol',      contopts.rel_tol,    ...
    'JPattern',    feval(handles{3},0,x0,parameters{:}) ...
  );



  t_end = 10;
  dydt = handles{2};
  f =@(t, y) dydt(t, y, parameters{:});
  disp("with dense Jac")
  tic
  sol = ode15s(f, [0 t_end], x0, int_opt);
  toc
  disp(sol.stats)
 
  dense_jac =  feval(handles{3},0,x0,parameters{:});
  


 
  disp("with Jpattern")
  tic
  sol = ode15s(f, [0 t_end], x0, int_opt);
  toc
  disp(sol.stats)
  
  

 

  odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d_sparse_jac', N));
  handles = feval(odefile);
  
  int_opt = odeset( ...
    'AbsTol',      contopts.abs_tol,    ...
    'RelTol',      contopts.rel_tol,    ...
    'Jacobian',    feval(handles{3},0,x0,parameters{:})' ...
  );

  disp("with sparse jac")
  tic
  sol = ode15s(f, [0 t_end], x0, int_opt);
  toc
  disp(sol.stats)
   
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

