function dx=f(t,x)
% right side
dx=[x(1)-x(2); x(2)-x(3); t*x(3)];
