% In this script we demonstrate how to find a stable equilibrium by time
% integration and continue it. The continuation is backwards with respect to the
% parameter L. The equilibrium becomes unstable after a branching point. The
% branching point occurs at L=pi.


handles = scalar_reaction_diffusion;
L = 5 * pi;
f = @(t,y) feval(handles{2}, t, y, L);

[t,y] = ode15s(f, [0 100], 0.2 * ones(1000,1));

figure
title('steps toward the stable equilibrium')
hold on;
for i = 1 : floor(size(y,1) / 10) : size(y,1)
  plot(y(i,:))
end

drawnow
equilibrium = y(end,:);
equilibrium = equilibrium(:);

[x0,v0] = init_EP_EP_L(@scalar_reaction_diffusion, equilibrium, L, 1);

options = contset( ...
  'MaxNumPoints',  18, ...
  'InitStepsize', 3, ...
  'MaxStepsize',  3, ...
  'Backward',     true, ...
  'Singularities', true, ...
  'FunTolerance',  1e-8, ...
  'VarTolerance',  1e-7, ...
  'MaxNewtonIters',       30, ...
  'max_rel_funcnorm_increase', 1e6, ...
  'MaxTestIters', 50, ...
  'contL_SmoothingAngle', 3, ...
  'contL_Loc_FunTolerance', 1e-7, ...
  'contL_Loc_VarTolerance', 1e-7 ...
);

figure
title('continuation steps')
hold on

s = contL(@equilibriumL, x0, v0, options, 'callback', @plot_equilibrium);


function plot_equilibrium(point, ~)
  plot(point.x(1:end-1)); 
  drawnow
end
