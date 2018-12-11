function bratu ( )
%
%    The Bratu equation includes a parameter lambda.  Depending on the
%    value of lambda, there may be 2, 1, or no solutions to the BVP.
%
%  Modified:   20 August 2017 by Mark Pekker
%
%    Original MATLAB version by Shampine, Kierzenka, Reichelt, John Burkardt (2013).
%
%  Reference:
%    Lawrence Shampine, Jacek Kierzenka, Mark Reichelt,
%    Solving boundary value problems for ordinary differential equations
%    in MATLAB with bvp4c.
%
  global lambda

  lambda_test = [ 0.45, 1.00, 3.50 ];

  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'BRATU:\n' );
  fprintf ( 1, '  MATLAB version\n' );
  fprintf ( 1, '  Use BVP4C to solve the following boundary value problem:\n' );
  fprintf ( 1, '  y" + lambda * exp ( y ) = 0\n' );
  fprintf ( 1, '  y(0) = 0, y(1) = 0\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  When lambda = 1, there are two solutions.\n' );
  fprintf ( 1, '  Try\n' );
  fprintf ( 1, '    y(x) = 0.1, y''(x) = 0.\n' );
  fprintf ( 1, '  and\n' );
  fprintf ( 1, '    y(x) = 3.0, y''(x) = 0.0\n' );

  figure_num = 0;
  
  t1 = tic;  
  for test = 1 : 3
    lambda = lambda_test(test);
    figure_num = bratu_solver ( lambda, figure_num );
  end 
  toc(t1)    
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'BRATU:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );
   
  return
end  

function figure_num = bratu_solver ( lambda, figure_num )

%    Input, real LAMBDA, the value of the parameter lambda.
%
%    Input/output, integer FIGURE_NUM, the number of figures being displayed.
%
%  Compute SOL1 with the first initial guess.
%
  %%x_init = linspace ( 0.0, 1.0, 5 );
  x_init = linspace ( 0.0, 1.0, 100 );
  y_init = [ 0.1, 0.0 ];
  solinit = bvpinit ( x_init, y_init );
 %% sol1 = bvp4c ( @bratu_ode, @bratu_bc, solinit );
   solver = 'bvp5c';
   %%solver = 'bvp4c';
%%end
  bvpsolver = fcnchk(solver);
  sol1 = bvpsolver ( @bratu_ode, @bratu_bc, solinit );
%
%  Compute SOL2 with the second initial guess.
%
  %%x_init = linspace ( 0.0, 1.0, 5 );
  x_init = linspace ( 0.0, 1.0, 100 );
  y_init = [ 3.0, 0.0 ];
  solinit = bvpinit ( x_init, y_init );
  %%sol2 = bvp4c ( @bratu_ode, @bratu_bc, solinit );
  sol2 = bvpsolver ( @bratu_ode, @bratu_bc, solinit );

  x = linspace ( 0.0, 1.0, 101 );
  y1 = deval ( sol1, x );
  y2 = deval ( sol2, x );
%
%  Display a plot of the two solutions.
%
  figure_num = figure_num + 1;
  figure ( figure_num )
  plot ( x, y1(1,:), 'r-', ...
         x, y2(1,:), 'g-', 'Linewidth', 2 );
  xlabel ( '<--- X --->' )
  ylabel ( '<--- Y --->' )
  title ( sprintf ( 'Bratu''s equation for \\lambda = %g\n', lambda ) );
  grid on
  filename = sprintf ( 'bratu_%f.png', lambda );
  print ( '-dpng', filename );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Saving plot file as "%s"\n', filename );

  return
end
function dydx = bratu_ode ( x, y )

%*****************************************************************************80
%
%% BRATU_ODE evaluates the right hand side of the ODE.
%
%    Input, real X, the point at which the ODE is to be evaluated.
%    Input, real Y(M), the value of the solution at X.
%    Output, real DYDX(M), the value of the right hand side given X and Y.
%
  global lambda

  dydx(1) = y(2);
  dydx(2) = - lambda * exp ( y(1) );

  return
end
function bc = bratu_bc ( ya, yb )
%
%% SAMPLE1_BC evaluates the boundary conditions.

%  Parameters:
%
%    Input, real YA(M), YB(M), the solution value at the left and right endpoints.
%    Output, real BC(2), the value of the boundary conditions.
%
  bc(1) = ya(1);
  bc(2) = yb(1);

  return
end
function timestamp ( )
%
%% TIMESTAMP prints the current YMDHMS date as a timestamp.
%
  t = now;
  c = datevec ( t );
  s = datestr ( c, 0 );
  fprintf ( 1, '%s\n', s );

  return
end
