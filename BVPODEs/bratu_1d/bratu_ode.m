function dydx = bratu_ode (x, y, lambda)
%
%    Input, real X, the point at which the ODE is to be evaluated.
%    Input, real Y(M), the value of the solution at X.
%    Output, real DYDX(M), the value of the right hand side given X and Y.
%
  dydx(1) = y(2);
  dydx(2) = - lambda * exp ( y(1) );
%---------------------------------------------------------