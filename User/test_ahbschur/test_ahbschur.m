function test_ahbschur
  % continuation of cycles cycles in fusion system
  clear global cds
  clear global lds
  N = 25;
  nphase = 3*(N-1);
  odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
  a = -1;
  b = -0.3;
  load('/home/carel/Documents/cl_matcontL/User/fusion_newton_picard/Data/fusion_np_from_previous_20-Mar-2019_18_43_33/point 3.mat', 'point')
  phi_0 = point.x(1:end-2);
  period = point.x(end-1);
  q_inf = point.x(end);
  
  global contopts cds
  contopts.abs_tol = 1e-12;
  contopts.rel_tol = 1e-12;
  
  parameters = {a ; b; q_inf};
  
  
  
  handles = feval(odefile);
  
  cds.dydt_ode = handles{2};
  cds.jacobian_ode = handles{3};
  cds.mv_count = 0;
  cds.integrator = @ode15s;

  int_opt = odeset( ...
    'AbsTol',      contopts.abs_tol,    ...
    'RelTol',      contopts.rel_tol,    ...
    'JPattern',    feval(handles{3},0,phi_0,parameters{:}) ...
  );

  
  
  cds.cycle_orbit = ode15s(...
    @(t, y) cds.dydt_ode(t, y, parameters{:}), ...
    [0 period], ...
    phi_0, ...
    int_opt);
  
  

  M = @(x) NewtonPicard.SingleShooting.monodromy_map(x, period, parameters, ...
                    contopts.abs_tol, contopts.rel_tol);
                  
  n_evals = 5;
  if false || true            
    tic
    [V,D] = eigs(M,nphase,n_evals);
    disp(diag(D))
    toc
    disp(cds.mv_count)
    cds.mv_count = 0;
  end
  
  
  options.k = n_evals;
  if false || false
    tic
    [q,t,flag] = ahbschur( @ahb_monodromy_map, ...
                      nphase, options, period, parameters, ...
                      contopts.abs_tol, contopts.rel_tol);
    disp(diag(t))
    toc
    disp(cds.mv_count)
  end
  options.BLSZ = 2;
  tic
  [q,t,flag] = ahbschur( @ahb_monodromy_map_parfor, ...
                    nphase, options, period, parameters, ...
                    cds.cycle_orbit, cds.jacobian_ode, ...
                    contopts.abs_tol, contopts.rel_tol);
  disp(diag(t))
  toc
  disp(cds.mv_count)
end

function result = ahb_monodromy_map(x, period, parameters, ...
                                                            abs_tol, rel_tol)
  result = zeros(size(x));
  for i=1:size(x,2)
    result(:,i) = NewtonPicard.SingleShooting.monodromy_map( ...
                                 x(:,i), period, parameters, abs_tol, rel_tol);
  end
end

function result = ahb_monodromy_map_parfor(x, period, parameters, ...
                                     cycle_orbit, jacobian_ode, abs_tol, rel_tol)
  result = zeros(size(x));
  global cds;
  cds.mv_count = cds.mv_count + size(x,2);
  parfor i=1:size(x,2)
    integration_opt = odeset('AbsTol', abs_tol, 'RelTol', rel_tol, ...
        'Jacobian',     @(t,y) feval(jacobian_ode, ...
                                t, deval(cycle_orbit,t), parameters{:}));
    dydt_mon = @(t, y) ...
      jacobian_ode(t, deval(cycle_orbit,t), parameters{:}) * y;
    [~,orbit] = ode15s(dydt_mon, [0 period], x(:,i), integration_opt);
    result(:,i) = orbit(end,:)';
  end
end

