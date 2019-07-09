function dfdx = jacobian(~,x,N,L,~,B,Dx,Dy) % unused inputs are t and A
  y  = x(N+1:2*N);
  x(N+1:2*N) =[];
  L2 = L^2;
  h  = 1/(N+1);
  cx = (Dx/L2)/(h*h);
  cy = (Dy/L2)/(h*h);
  % Sparse jacobian
  A = zeros(2*N,3);
  A(1:N-1,2)=cx;
  A(1:N,3)=-2*cx -(B+1) + 2*x(1:N).*y(1:N);
  A(1:N,4)=cx;
  A(N+1:2*N,2) = cy;
  A(N+1:2*N,3) = -2*cy -x(:).*x(:);
  A(N+2:2*N,4) = cy;
  A(1:N,1) = B - 2*x(:).*y(:);
  A(N+1:2*N,5) = x(:).*x(:);
  dfdx = spdiags(A, [-N,-1:1,N] , 2*N, 2*N);