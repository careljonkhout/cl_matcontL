function out = bratu_1d
%
% Odefile of 1-d Brusselator model ?
% This file uses a non-constant approximate solution as a starting guess. ?

out{1} = @init;
out{2} = @fun_eval;
out{3} = [];@jacobian;
out{4} = [];%@jacobianp;
out{5} = [];@usernorm;
out{6} = [];@user1;
% ----------------------------------------------------------------------
function dxdt = fun_eval(t,x,N,lambda)
y  = x(N+1:2*N);
x(N+1:2*N) =[];
x0 = A; x1 = A;
y0 = B/A; y1 = B/A;
L2 = L^2;
h  = 1/(N+1);
cx = (Dx/L2)/(h*h);
cy = (Dy/L2)/(h*h);
dxdt = zeros(N,1);
dydt = zeros(N,1);
dxdt(1) = (x0-2*x(1)+x(2))*cx + A - (B+1)*x(1) + x(1)*x(1)*y(1);
dxdt(N) = (x(N-1)-2*x(N)+x1)*cx + A - (B+1)*x(N) + x(N)*x(N)*y(N);
dydt(1) = (y0-2*y(1)+y(2))*cy + B*x(1) - x(1)*x(1)*y(1);
dydt(N) = (y(N-1)-2*y(N)+y1)*cy + B*x(N) - x(N)*x(N)*y(N);
i=2:N-1;
dxdt(i) = (x(i-1)-2*x(i)+x(i+1))*cx + A - (B+1)*x(i) + x(i).*x(i).*y(i);
dydt(i) = (y(i-1)-2*y(i)+y(i+1))*cy + B*x(i) - x(i).*x(i).*y(i);
dxdt = [dxdt; dydt];
% --------------------------------------------------------------------------
function y0 = init(N,lambda)

x_init = linspace (0.0, 1.0, N);
y_init = [0.1, 0.0];
solinit = bvpinit (x_init, y_init);
y0_1 = solinit.y(1,1:N)';
y0_2 = solinit.y(2,1:N)';

y0 = zeros(2*N,1);
y0 = [y0_1; y0_2];
% --------------------------------------------------------------------------
function dfdx = jacobian(t,x,N,lambda)
y  = x(N+1:2*N);
x(N+1:2*N) =[];
x0 = A; x1 = A;
y0 = B/A; y1 = B/A;
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
% --------------------------------------------------------------------------
function dfdp = jacobianp(t,x,N,lambda)
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
% -------------------------------------------------------------------------
function normuser = usernorm(arg)
n = length(arg);
normuser = sqrt(1/n)*norm(arg);
% -------------------------------------------------------------------------
function uf1 = user1(t,x,N,lambda)
uf1 = L - 6.5e-002;
% ---------------------------------------------------------------------
% ----