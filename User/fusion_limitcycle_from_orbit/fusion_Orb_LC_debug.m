load('fusion_2_cycles_ready_for_initOrbLC')

p1 = [a; b; -0.72];
p2 = [a; b; -0.72 + 0.001];
 
[x0_2,~] = initOrbLC(odefile, t2, x2, p1, ap, ntst, ncol, tolerance);
[x0_1,~] = initOrbLC(odefile, t3, x3, p2, ap, ntst, ncol, tolerance);
v0 = x0_2 - x0_1;

[sout, datafile] = contL(@limitcycleL,x0_1,v0,opt); %#ok<*ASGLU>

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

