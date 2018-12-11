function single_shooting_fusion
  % continuation of cycles cycles in fusion
  % slower than orthogonal collocation
  
  N = 25;
  nPhases = 3*(N-1);
  odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
  a = -1;
  b = -0.3;
  q_inf = -0.72;
  parameters = {a;b;q_inf};
  poincare_tolerance = 1e-4;
  
  
  handles = feval(odefile);
  
  int_opt = odeset( ...
    'AbsTol',      1e-10,    ...
    'RelTol',      1e-13,    ...
    'Jacobian',     @(t,y) feval(handles{3},t,y,parameters{:}) ...
  );


  x0 = ones(nPhases,1);
  dydt = handles{2};
  f =@(t, y) dydt(t, y, parameters{:});
  [t1, x1] = ode15s(f, [0 150], x0, int_opt);

  draw_plots = false || false;
  if draw_plots
    figure(1)
    plot(t1,x1)
    %title(sprintf());
    %xlabel('t');
    %ylabel('todo: add label');
  end

  approximate_period = 12;
  
  x0 = x1(end,:)';
  v0 = f(0,x0)';
  int_opt = odeset(int_opt, 'Events', @returnToPlane);
  global cds;
  cds.nDiscretizationPoints = 100;
  [t2,x2] = ode15s(f, 0:approximate_period, x1(end,:), int_opt); 
  T = t2(end);
  int_opt = odeset(int_opt, 'Events', []);
  [t3,x3] = ode15s(f, linspace(0,T,cds.nDiscretizationPoints), x2(end,:), int_opt); 
  
  cds.x0_prime = zeros(cds.nDiscretizationPoints,nPhases);
  for i=1:cds.nDiscretizationPoints
    cds.x0_prime(i,:) = f(0,x3(i,:));
  end
  


  % use poincare section to converge to 
  %[t3,x3] = ode15s(f, 0:5*approximate_period,x1(end,:)',integration_opt);


  if draw_plots
    figure(2)
    hold on;
    plot(t2,x2)
     title(sprintf( ...
      'brusselator N:2 L:%.2f a:%.2f n:%.2f Dx:%.2f Dy:.2f', ...
      parameters{:}));
    xlabel('t')
    ylabel('x_1,x_2,y_1,y_2');
  end

  if draw_plots
    figure(3)
    hold on;
    plot(t3,x3)
    title(sprintf( ...
      'brusselator N:2 L:%.2f a:%.2f n:%.2f Dx:%.2f Dy:.2f', ...
      parameters{:}));
    xlabel('t')
    ylabel('n_1, ..., n_{N-1},U_1, ..., U_{N-1},z_1,...,z_{N-1}');
  end

  
  if draw_plots
    figure(3)
    plot(x2(:,1),x2(:,2))
   title(sprintf( ...
      'brusselator N:2 L:%.2f a:%.2f n:%.2f Dx:%.2f Dy:.2f', ...
      parameters{:}));
    xlabel('n_1')
    ylabel('U_1')

    drawnow
  end

  %% Continue limit cycle from orbit

  

  opt = contset();
  opt = contset(opt, 'MaxNumPoints',   3);
  opt = contset(opt, 'InitStepsize',   1e-2);
  opt = contset(opt, 'MinStepsize',    1e-4);
  opt = contset(opt, 'MaxStepsize',    1e-1);
  opt = contset(opt, 'MaxNewtonIters', 3);
  opt = contset(opt, 'MaxCorrIters',   10);
  opt = contset(opt, 'MaxTestIters',   10);
  opt = contset(opt, 'VarTolerance',   1e-6);
  opt = contset(opt, 'FunTolerance',   1e-6);
  % we don't want to adapt
  % since it is not implemented
  opt = contset(opt, 'Adapt',          1000*1000*1000);
  opt = contset(opt, 'MaxNumPoints',   100);
  opt = contset(opt, 'CheckClosed',    50);
  opt = contset(opt, 'Multipliers',    true);
  opt = contset(opt, 'Backward',       false);
  opt = contset(opt, 'Singularities',  false);
  opt = contset(opt, 'CIS_UsingCIS',   false);

   
  cds.probfile = odefile;
  cds.nap = 1;
  cds.ActiveParams = 3;
  cds.ndim = nPhases + cds.nap + 1;
  cds.P0 = cell2mat(parameters);
  cds.options = contset();
  cds.symjac = false;
  cds.usernorm = [];
  cds.probfile = odefile;
  [s, datafile] = contL(@single_shooting,[x0; T; cds.P0(cds.ActiveParams)],[],opt); 


  figure
  hold on;
  coordinate1 = 1;
  coordinate2 = 2;
   title(sprintf( ...
      'brusselator N:%d a:%.2f n:%.2f q_inf:%.2f', ...
       N,a,b,q_inf));
   xlabel('x_1')
   ylabel('y_1')
  [x, v, h, mult] = loadPoint(datafile); % DV: load computed cycles
  %load('Data\testbruss_Orb_LC.mat')    % DV: load singular points

  global lds
  close all
  figure
  hold on
  for i=1:10:opt.MaxNumPoints
    xx = x(1:end-1-length(lds.ActiveParams),i);
    coord1_vals = xx(coordinate1:lds.nphase:end);
    coord2_vals = xx(coordinate2:lds.nphase:end);
    plot(coord1_vals,coord2_vals,'b');
    title(sprintf( ...
      'brusselator N:%d a:%.2f n:%.2f q_{inf}:varied', ...
       N,parameters{:}));
    xlabel('n_1')
    ylabel('U_1')
  end
  figure
  plot(x(end,:),x(end-1,:))
  title(sprintf( ...
      'brusselator N:%d a:%.2f n:%.2f q_{inf}:varied', ...
       N,parameters{:}));
  xlabel('q_{inf}')
  ylabel('period')

  figure
  hold on;
  nMults = size(mult,1);
  for i=nMults-10:nMults
    plot(x(end,:),mult(i, :))
  end
  
  function [value, isterminal, direction] = returnToPlane(t, x)
     % v0 should be a row vector
     % x and x0 should be column vectors
     value = v0*(x-x0);
     isterminal =  sum((x-x0).^2) < poincare_tolerance && t > 1;
     direction = 1;
  end

end
