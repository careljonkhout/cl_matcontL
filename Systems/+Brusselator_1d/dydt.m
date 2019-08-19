function dydt = dydt(~,x_and_y,N,L,A,B,Dx,Dy) % unused input is t
  x  = x_and_y(1:N);    % 
  y  = x_and_y(N+1:end);
  x0 = A; x1 = A;
  y0 = B/A; y1 = B/A;
  L2 = L^2;
  h  = 1/(N+1);
  cx = (Dx/L2)/(h*h);
  cy = (Dy/L2)/(h*h);
  dxdt(1) = (x0     - 2*x(1) + x(2))*cx + A - (B+1)*x(1) + x(1)*x(1)*y(1);
  dydt(1) = (y0     - 2*y(1) + y(2))*cy + B*x(1) - x(1)*x(1)*y(1);
  dxdt(N) = (x(N-1) - 2*x(N) + x1  )*cx + A - (B+1)*x(N) + x(N)*x(N)*y(N);
  dydt(N) = (y(N-1) - 2*y(N) + y1  )*cy + B*x(N) - x(N)*x(N)*y(N);
  
  i=2:N-1;
  dxdt(i) = (x(i-1)-2*x(i)+x(i+1))*cx + A - (B+1)*x(i) + x(i).*x(i).*y(i);
  dydt(i) = (y(i-1)-2*y(i)+y(i+1))*cy + B*x(i) - x(i).*x(i).*y(i);
  dydt = [dxdt(:); dydt(:)];