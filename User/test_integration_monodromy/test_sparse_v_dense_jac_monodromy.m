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
  tol = 1e-13;
  contopts.abs_tol = tol;
  contopts.rel_tol = tol;
  contopts.monodromy_map_abs_tol = tol;
  contopts.monodromy_map_rel_tol = tol;
  contopts.integration_abs_tol   = tol;
  contopts.integration_rel_tol   = tol;
  cds.nphases = 3*(N-1);

  handles = feval(odefile);

  f =@(t, y) feval(handles{2},t, y, parameters{:});
  
  cds.integrator = @ode15s;
  cds.jacobian_ode = handles{3};
  
  pattern = feval(handles{3},0,ones(cds.nphases,1),parameters{:});
    integration_opt = odeset(...
    'AbsTol',      contopts.abs_tol,    ...
    'RelTol',      contopts.rel_tol    ...
  );

  
  cds.cycle_orbit  = ode15s(f, [0 period], x0, integration_opt);
  cds.mv_count = 0;
  if false 
    monodromy_map = @(x) monodromy_map_finite_differences( ...
      x0, x, period, parameters, 1e-4);
    disp('by finite differences')
    tic
    [~, multiplier_matrix] = eigs(monodromy_map, cds.nphases,5);
    toc
    multipliers = diag(multiplier_matrix);
    disp(multipliers);
  end
  integration_opt = odeset(...
    'AbsTol',      contopts.abs_tol,    ...
    'RelTol',      contopts.rel_tol,    ...
    'Jacobian',     @(t,y) feval(handles{3},t,y,parameters{:}) ...
  );
  
  monodromy_map = @(x) NewtonPicard.SingleShooting.monodromy_map( ...
    x, period, parameters);
  
  disp('by integrating dense jacobian');
  tic
  Mx = monodromy_map(f(0,x0));
  toc
  disp(norm(Mx - f(0,x0)));


 
  odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d_sparse_jac', N));
  handles = feval(odefile);
  
  
  monodromy_map = @(x) monodromy_map_from_jac(x, period, parameters);
  
  
  disp('by integrating sparse jacobian');
  tic
  Mx = monodromy_map(f(0,x0));
  toc
  disp(norm(Mx - f(0,x0)));

  pattern = feval(handles{3},0,ones(3*(N-1)),parameters{:})';
  monodromy_map = @(x) monodromy_map_pattern(x, period, parameters);
  disp('with JPattern');
  tic
  Mx = monodromy_map(f(0,x0));
  toc
  disp(norm(Mx - f(0,x0)));
 % figure
%  plot(Mx ./ f(0,x0));
  
  integration_opt = odeset(...
    'AbsTol',      contopts.abs_tol/2,    ...
    'RelTol',      contopts.rel_tol/2   ...
  );
  monodromy_map = @(x) monodromy_map_finite_differences(x0, x, period, parameters,1e-4);
  disp('with finite differences');
  tic
  Mx = monodromy_map(f(0,x0));
  toc
  disp(norm(Mx - f(0,x0)));
  
  integration_opt = odeset(...
    'AbsTol',      contopts.abs_tol/2,    ...
    'RelTol',      contopts.rel_tol/2   ...
  );

  monodromy_map = @(x) monodromy_map_finite_differences_4(x0, x, period, parameters,1e-3);
  disp('with finite differences 4 h=1e-3');
  tic
  Mx = monodromy_map(f(0,x0));
  toc
  disp(norm(Mx - f(0,x0)));

  monodromy_map = @(x) monodromy_map_finite_differences_4(x0, x, period, parameters,1e-4);
  disp('with finite differences 4 h=1e-4');
  tic
  Mx = monodromy_map(f(0,x0));
  toc
  disp(norm(Mx - f(0,x0)));

 
  monodromy_map = @(x) monodromy_map_finite_differences_4(x0, x, period, parameters,1e-5);
  disp('with finite differences 4 h=1e-5');
  tic
  Mx = monodromy_map(f(0,x0));
  toc
  disp(norm(Mx - f(0,x0)));

  

  
%   function x_end = shoot(x, period, ~)
%     [~, orbit] = ode15s(f, [0 period], x, integration_opt);
%     x_end = orbit(end,:)';
%   end



  function d_phi__d_x = monodromy_map_finite_differences(...
      x_cycle, x, period, parameters, stepsize)
    h = stepsize;
    f = @(t, x) feval(handles{2},0,x,parameters{:});
    ff = @(t, x1_and_x2) [f(0, x1_and_x2(1:cds.nphases));
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

  function d_phi__d_x = monodromy_map_finite_differences_4(...
      x_cycle, x, period, parameters, stepsize)
    
    h = stepsize;
    f = @(t, x) feval(handles{2},0,x,parameters{:});
    ff = @(t, x1234) [f(0, x1234(1                :    cds.nphases))
                      f(0, x1234(    cds.nphases+1:2 * cds.nphases))
                      f(0, x1234(2 * cds.nphases+1:3 * cds.nphases))
                      f(0, x1234(3 * cds.nphases+1:end))];
    [~, orbit] = ode15s(ff, ...
      [0 period], ...
      [x_cycle - 2 * h * x; x_cycle -     h * x; ...
       x_cycle +     h * x; x_cycle + 2 * h * x], ...
      integration_opt);
    phi_x1234 = orbit(end,:)';
    phi_x1 = phi_x1234(              1:  cds.nphases);
    phi_x2 = phi_x1234(  cds.nphases+1:2*cds.nphases);
    phi_x3 = phi_x1234(2*cds.nphases+1:3*cds.nphases);
    phi_x4 = phi_x1234(3*cds.nphases+1:end          );
    
    d_phi__d_x = (phi_x1 - 8 * phi_x2 + 8 * phi_x3 - phi_x4 )/h/12; 
  end

  function Mx = monodromy_map_from_jac(phases_0, period, parameters)
    int_opt = odeset(...
      'AbsTol',       1e-13, ...
      'RelTol',       1e-13, ...
      'Jacobian',     @(t,y) feval(handles{3}, ...
                        t, deval(cds.cycle_orbit,t), parameters{:})' ...
    );                
    dydt_mon = @(t, y) ...
      feval(handles{3},t, deval(cds.cycle_orbit,t), parameters{:})' * y;
    [~,orbit] = ode15s(dydt_mon, [0 period], phases_0, int_opt);
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

