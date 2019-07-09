% note this Jacobian for parameters is incomplete
function dfdp = jacobian_params(~,x,N,L,A,B,Dx,Dy) % unused parameter is t
  y  = x(N+1:2*N);
  x(N+1:2*N) =[];
  x0 = A; x1 = A;
  y0 = B/A; y1 = B/A;
  L2 = L^2;
  h = 1/(N+1);
  cx = (Dx/L2)/(h*h);
  cy = (Dy/L2)/(h*h);
  kx = (-2/L)*cx;
  ky = (-2/L)*cy;
  Sx = zeros(N,1);
  Sy = zeros(N,1);
  Sx(1) = kx*(x0-2*x(1)+x(2));
  Sy(1) = ky*(y0-2*y(1)+y(2));
  Sx(N) = kx*(x(N-1)-2*x(N)+x1);
  Sy(N) = ky*(y(N-1)-2*y(N)+y1);
  i=2:N-1;
  Sx(i) = kx*(x(i-1)-2*x(i)+x(i+1));
  Sy(i) = ky*(y(i-1)-2*y(i)+y(i+1));
  dfdp = [ zeros(2*N,1) [Sx;Sy] ];