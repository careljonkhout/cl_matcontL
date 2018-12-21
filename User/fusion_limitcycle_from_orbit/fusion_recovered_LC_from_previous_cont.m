run_init_if_needed
clc
if ~ exist('x75','var')
  [x75,v75,~,mult75] = loadPoint(...
    fullfile('Data','fusion_Orb_LC_N_75_16-Dec-2018_16_58_59.dat'));
end
nPhases = 3*(N-1);
x0 = x75(:,end);
v0 = v75(:,end);




N = 75;                     
odefile = str2func(sprintf('fusion_precomputed_with_sage_N_%d', N));
a = -1;
b = -0.3;


period = x0(end-1);
q_inf = x0(end);
parameters = {a ; b; q_inf};


% continuation of limit cycles in fusion system
global poincare_tolerance
poincare_tolerance = 1e-4;

handles = feval(odefile);

integration_opt = odeset( ...
  'AbsTol',      1e-13,    ...
  'RelTol',      1e-13,    ...
  'Jacobian',    @(t,y) feval(handles{3},t,y,parameters{:}), ...
  'Events',      @returnToPlane ...
);


global plane

x = x0(1:nPhases);
for i=1:(ntst*ncol)
  plane.x0 = x0((i+1)*nPhases+1:(i+2)*nPhases);
  plane.v0 = v0((i+1)*nPhases+1:(i+2)*nPhases);
  [t, trajectory] = ode15s(f, linspace(0,period,1000), x, integration_opt);
  x = trajectory(end,:)';
  fprintf('%.5f\n',t);
end
return

dydt = handles{2};
f =@(t, y) dydt(t, y, parameters{:});


draw_plots = false || true;
if draw_plots
  figure(1)
  plot(t1,x1)
  title(sprintf( ...
    'fusion N:%d a:%.2f n:%.2f q_{inf}:%.2f', ...
    N,a,b,q_inf));
  xlabel('t');
  ylabel('n_1, ..., n_{N-1},U_1, ..., U_{N-1},z_1,...,z_{N-1}');
end


if draw_plots
  figure(3)
  plot(x1(:,1),x1(:,2))
  title(sprintf( ...
    'fusion N:%d a:%.2f n:%.2f q_{inf}:%.2f', ...
     N,a,b,q_inf));
  xlabel('n_1')
  ylabel('U_1')

  drawnow
end

%% Continue limit cycle from orbit
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
opt = contset(opt, 'Adapt',          3);
opt = contset(opt, 'MaxNumPoints',   100);
opt = contset(opt, 'CheckClosed',    50);
opt = contset(opt, 'Multipliers',    true);
opt = contset(opt, 'Backward',       false);
opt = contset(opt, 'Singularities',  false);
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

 function [value, isterminal, direction] = returnToPlane(t, x)
   % v0, x, and x0 should be column vectors
   global plane
   global poincare_tolerance
   value = plane.v0'*(x-plane.x0);
   isterminal =  sum((x-plane.x0).^2) < poincare_tolerance && t > 1;
   direction = 1;
 end
