close all

f = @(t, x) [x(2); -x(1)];
t = linspace(0, 2 * pi, 2);

sol_many_points = ode15s(f, t, [0 1]);
plot(sol.x, sol.y);