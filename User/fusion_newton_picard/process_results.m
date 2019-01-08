datafile=fullfile('Data','fusion_newton_picard_03-Jan-2019_18_04_18.dat');
 x = loadPoint(datafile); % DV: load computed cycles
plotcycles = false || true;

N=25;
a=1.0;
b=-0.3;
if plotcycles

  close all
  figure
  hold on
  for i=203:203
    xx = x(1:end,i);
    period                       = xx(end-1);
    phases_0                     = xx(1:end-2);
    parameters                   = cds.P0;
    parameters(cds.ActiveParams) = xx(end);
    parameters                   = num2cell(parameters);
    [t,y] = compute_cycle(phases_0,period,parameters);
    coord1_vals = y(:,1);
    coord2_vals = y(:,3);
    plot(coord1_vals,coord2_vals,'b');
    title(sprintf('fusion N:%d a:%.2f b:%.2f',N,a,b));
    xlabel('n_1')
    ylabel('U_1')
  end
end

figure

title(sprintf('fusion N:%d a:%.2f b:%.2f',N,a,b));
xlabel('q_{inf}')
ylabel('period')
plot(x(end,:),x(end-1,:));
figure
plot(x(end,:));
figure
plot(x(end-1,:));

function [value, isterminal, direction] = returnToPlane(t, x)
  global cds;
  % x and should be a column vectors
  value = cds.dydt_0'*(x-cds.x0);
  isterminal = t > 1 && sum((x-cds.x0).^2) < cds.poincare_tolerance;
  direction = 1;
end

function [t,x] = compute_cycle(x, period, parameters)
  global cds
  f =@(t, y) cds.dydt_ode(t, y, parameters{:});
  integration_opt = odeset(...
    'AbsTol',      1e-10,    ...
    'RelTol',      1e-10,    ...
    'BDF',         'off',   ...
    'MaxOrder',     5,      ...
    'NormControl',  'off',  ...
    'Refine',       1,      ...
    'Jacobian',     @(t,y) feval(cds.jacobian_ode,t,y,parameters{:}) ...
  );
  [t,x] = ode15s(f, [0 period], x, integration_opt);
end

