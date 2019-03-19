function out = cubic_quintic_Ginzburg_Landau
%
% Odefile of 1-d cubic quintic Ginzburg Landau equation
% discretized using a finite differences on a equidistant mesh

out{1} = @init;
out{2} = @fun_eval;
out{3} = [];%@jacobian;
out{4} = [];%@jacobianp;
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
out{10} = @usernorm;
out{11} = @user1;
% ----------------------------------------------------------------------
function du = fun_eval(~,x,r,c3,c5,mu,nu) %ignored parameter is t (time)
N = length(x) / 2;
h = 2 * pi / (N-1);
u = reshape(x,2,N);
du = zeros(size(u));
for i=1:N
  du(:,i) = reaction_term(u(:,i),r,c3,c5,mu,nu);
end
% discretized zero Neumann BC left:
du(:,1) = du(:,1) + (- 2 * u(:,1) + 2 * u(:, 2  )) / h^2;
% discretized zero Neumann BC right:
du(:,N) = du(:,N) + (- 2 * u(:,N) + 2 * u(:, N-1)) / h^2;

for i=2:N-1
  du(:,i) = du(:,i) + (u(:,i-1) - 2 * du(:,i) + u(:,i+1))/(h^2);
end

function f = reaction_term(u,r,c3,c5,mu,nu)
  modulus_squared = u(1)^2 + u(2)^2;
  f = ([ r -nu ; nu r ] - modulus_squared * [ c3 - mu ; mu c3 ]) * u  ...
      - c5 * modulus_squared^2 *u;
% --------------------------------------------------------------------------
function init()


% --------------------------------------------------------------------------
function dfdx = jacobian(t,x,N,L,A,B,Dx,Dy)
% y  = x(N+1:2*N);
% x(N+1:2*N) =[];
% x0 = A; x1 = A;
% y0 = B/A; y1 = B/A;
% L2 = L^2;
% h  = 1/(N+1);
% cx = (Dx/L2)/(h*h);
% cy = (Dy/L2)/(h*h);
% % Sparse jacobian
% A = zeros(2*N,3);
% A(1:N-1,2)=cx;
% A(1:N,3)=-2*cx -(B+1) + 2*x(1:N).*y(1:N);
% A(1:N,4)=cx;
% A(N+1:2*N,2) = cy;
% A(N+1:2*N,3) = -2*cy -x(:).*x(:);
% A(N+2:2*N,4) = cy;
% A(1:N,1) = B - 2*x(:).*y(:);
% A(N+1:2*N,5) = x(:).*x(:);
% dfdx = spdiags(A, [-N,-1:1,N] , 2*N, 2*N);
% --------------------------------------------------------------------------
function dfdp = jacobianp(t,x,N,L,A,B,Dx,Dy)
% y  = x(N+1:2*N);
% x(N+1:2*N) =[];
% x0 = A; x1 = A;
% y0 = B/A; y1 = B/A;
% L2 = L^2;
% h = 1/(N+1);
% cx = (Dx/L2)/(h*h);
% cy = (Dy/L2)/(h*h);
% kx = (-2/L)*cx;
% ky = (-2/L)*cy;
% Sx = zeros(N,1);
% Sy = zeros(N,1);
% Sx(1) = kx*(x0-2*x(1)+x(2));
% Sy(1) = ky*(y0-2*y(1)+y(2));
% Sx(N) = kx*(x(N-1)-2*x(N)+x1);
% Sy(N) = ky*(y(N-1)-2*y(N)+y1);
% i=2:N-1;
% Sx(i) = kx*(x(i-1)-2*x(i)+x(i+1));
% Sy(i) = ky*(y(i-1)-2*y(i)+y(i+1));
% dfdp = [ zeros(2*N,1) [Sx;Sy] ];
% -------------------------------------------------------------------------

