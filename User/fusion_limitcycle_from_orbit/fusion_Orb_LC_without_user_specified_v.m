if ~ exist('contL', 'file')
  fullpath = mfilename('fullpath');
  file_directory = fullpath(1:end-length(mfilename));
  cd([file_directory '../..']);
  init
  cd(file_directory)
end



% continuation of cycles cycles in brusselator
N = 25;                     
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
a = -1;
b = -0.3;
q_inf = -0.72;
poincare_tolerance = 1e-4;
parameters = {a ; b; q_inf};
handles = feval(odefile);

integration_opt = odeset( ...
  'AbsTol',      1e-10,    ...
  'RelTol',      1e-13,    ...
  'Jacobian',     @(t,y) feval(handles{3},t,y,parameters{:}) ...
);


x0 = ones(3*(N-1),1);
dydt = handles{2};
f =@(t, y) dydt(t, y, parameters{:});
[t1, x1] = ode15s(f, [0 150], x0, integration_opt);

draw_plots = false || true;
if draw_plots
  figure(1)
  plot(t1,x1)
  title(sprintf( ...
    'fusion N:%d a:%.2f n:%.2f q_{inf}:%.2f', ...
    N,a,b,q_inf));
  xlabel('t');
  ylabel('n_1, ..., n_{N-1},U_1, ..., y_{N-1},z_1,...,z_{N-1}');
end

approximate_period = 10;

x0 = x1(end,:)';
v0 = f(0,x0)';
%  integration_opt = odeset(integration_opt, 'Events', @returnToPlane);

[t2,x2] = ode15s(f, 0:0.01:approximate_period, x1(end,:), ...
  integration_opt); 



if draw_plots && true
  figure(2)
  hold on;
  plot(t2,x2)
  title(sprintf( ...
    'fusion N:%d a:%.2f n:%.2f q_{inf}:%.2f', ...
     N,a,b,q_inf));
  xlabel('t')
  ylabel('n_1, ..., n_{N-1},U_1, ..., U_{N-1},z_1,...,z_{N-1}');
end



if draw_plots
  figure(3)
  plot(x2(:,1),x2(:,2))
  title(sprintf( ...
    'fusion N:%d a:%.2f n:%.2f q_{inf}:%.2f', ...
     N,a,b,q_inf));
  xlabel('n_1')
  ylabel('U_1')

  drawnow
end

%% Continue limit cycle from orbit
q_inf = -0.72;
parameters = {a ; b; q_inf};
p = cell2mat(parameters);
ntst = 20;
ncol = 4;
tolerance = 1e-2;

opt = contset();
opt = contset(opt, 'MaxNumPoints',   3);
opt = contset(opt, 'InitStepsize',   3);
opt = contset(opt, 'MinStepsize',    1e-4);
opt = contset(opt, 'MaxStepsize',    20);
opt = contset(opt, 'MaxNewtonIters', 3);
opt = contset(opt, 'MaxCorrIters',   10);
opt = contset(opt, 'MaxTestIters',   10);
opt = contset(opt, 'VarTolerance',   1e-6);
opt = contset(opt, 'FunTolerance',   1e-6);
opt = contset(opt, 'TestTolerance',  1e-5);
opt = contset(opt, 'Adapt',          3);
opt = contset(opt, 'MaxNumPoints',   100);
opt = contset(opt, 'CheckClosed',    50);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       false);
opt = contset(opt, 'Singularities',  true);
opt = contset(opt, 'CIS_UsingCIS',   false);

ap = 3;



[x0_2,~] = initOrbLC(odefile, t2, x2, p, ap, ntst, ncol, tolerance);

[sout, datafile] = contL(@limitcycleL,x0_2,[],opt); %#ok<*ASGLU>


figure
hold on;
coordinate1 = 1;
coordinate2 = 2;
 title(sprintf( ...
    'fusion N:%d a:%.2f n:%.2f q_inf:%.2f', ...
     N,a,b,q_inf));
 xlabel('n_1')
 ylabel('U_1')
[x, v, h, mult] = loadPoint(datafile); % DV: load computed cycles


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
    'fusion N:%d a:%.2f n:%.2f q_{inf}:varied', ...
     N,a,b));
  xlabel('n_1')
  ylabel('U_1')
end
figure
plot(x(end,:),x(end-1,:))
title(sprintf( ...
    'fusion N:%d a:%.2f n:%.2f q_{inf}:varied', ...
     N,a,b));
xlabel('q_{inf}')
ylabel('period')

figure
hold on;
nMults = size(mult,1);
for i=nMults-10:nMults
  plot(x(end,:),mult(i, :))
end

