function test_monodromy_finite_differences
  clear global cds
  clear global lds
  N = 25;                     
  odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
  handles = feval(odefile);
  a = -1;
  b = -0.3;
  q_inf = -0.72;
  parameters = {a ; b; q_inf};
  load('/home/carel/Documents/cl_matcontL/User/fusion_newton_picard/Data/fusion_np_from_previous_20-Mar-2019_18_43_33/point 3.mat', 'point')
  x0 = point.x(1:end-2);
  period = point.x(end-1);
  parameters{3} = point.x(end);
  global contopts cds
  tol = 1e-13;
  contopts.abs_tol = tol;
  contopts.rel_tol = tol;
  contopts.monodromy_map_abs_tol = tol;
  contopts.monodromy_map_rel_tol = tol;
  
  pattern = feval(handles{3},0,ones(3*(N-1)),parameters{:});
  pattern = pattern ~= 0;
  pattern = sparse(pattern);
  
  double_pattern = blkdiag(pattern, pattern);



  
  
  cds.nphases = 3*(N-1);

  handles = feval(odefile);

  f =@(t, y) feval(handles{2},t, y, parameters{:});
  
  cds.integrator = @ode15s;
  cds.jacobian_ode = handles{3};
 
 % cds.cycle_orbit  = ode15s(f, [0 period], x0, integration_opt);
  

  odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d_sparse_jac', N));
  handles = feval(odefile);
  

  global errors
  stepsizes = power(10,-2:-0.25:-10);
  errors = zeros(size(stepsizes));
  
  i = 1;
  
    integration_opt = odeset(...
      'AbsTol',       contopts.monodromy_map_abs_tol, ...
      'RelTol',       contopts.monodromy_map_rel_tol, ...
      'JPattern',     double_pattern ...
  ); 
  
  for stepsize = stepsizes
    monodromy_map = @(x) monodromy_map_finite_differences(x0, x, period, parameters, stepsize);
    tic
    Mx =  monodromy_map(f(0,x0));
    toc
    errors(i) = norm(Mx - f(0,x0));
    fprintf('%d %.2e %.4e\n', i, stepsize, errors(i));
    i = i + 1;
  end
  plot(stepsizes, errors);
  


  function d_phi__d_x = monodromy_map_finite_differences(...
      x_cycle, x, period, parameters, stepsize)
    h = stepsize;
    f = @(t, x) feval(handles{2},0,x,parameters{:});
    ff = @(t, x1_and_x2) [f(0, x1_and_x2(1:cds.nphases))
                          f(0, x1_and_x2(cds.nphases+1:end))];
    [~, orbit] = ode15s(ff, ...
      [0 period], ...
      [x_cycle - h * x; x_cycle+h * x], ...
      integration_opt);
    
    phi_x1__and__phi_x2 = orbit(end,:)';
    phi_x1 = phi_x1__and__phi_x2(1:cds.nphases);
    phi_x2 = phi_x1__and__phi_x2(cds.nphases+1:end);
    
    d_phi__d_x = (phi_x2 - phi_x1)/h/2; 
  end
end

