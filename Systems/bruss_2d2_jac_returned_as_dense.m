
function out = bruss_2d2_jac_returned_as_dense
%
% Odefile of 2-d Brusselator model, 
% C.S. Chien, Z.~Mei, and C.L. Shen.???
% Numerical continuation at double bifurcation points of a
% reaction-diffusion problem, Int. J. Bifur. and Chaos, 8(1):117--139, 1997.??
%

out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = [];%@jacobianp;
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
out{10} = @usernorm;
%%out{6} = [];%@userfun1;
out{11} = @user1;
out{12} = @user2;
out{13} = @user3;
out{14} = @user4;

% --------------------------------------------------------------------------

%%function dfdt = fun_eval(t,y,N,A,B)
function dfdt = fun_eval(t,y,N,L,A,B,Dx,Dy)

%%dfdt = fj(y,N,A,B,1);
dfdt = fj(y,N,L,A,B,Dx,Dy,1);

% --------------------------------------------------------------------------

function x0 = init(N,L,A,B,Dx,Dy)     % CIS
n = N*N;
x0 = zeros(2*n,1);

for i=1:n
  %%x0(i)   = A;
  %%x0(n+i) = B/A;
  y0(i) = 0;
  y0(n+i) = 0;
end


% -----------------------------------------------------------------------

%%function dfdxy = jacobian(t,y,N,A,B)
function dfdxy = jacobian(t,y,N,L,A,B,Dx,Dy)

%%dfdxy = fj(y,N,A,B,2);    % Sparse jacobian
dfdxy = full(fj(y,N,L,A,B,Dx,Dy,2));    % Sparse jacobian

% -----------------------------------------------------------------------

%%function dfdp = jacobianp(t,y,N,A,B)
function dfdp = jacobianp(t,x,N,L,A,B,Dx,Dy)
% df/dB

x = y(1:N);


Sx = zeros(N,1);
Sy = zeros(N,1);

Sx(1) = -x(1);
Sx(N) = -x(N);

Sy(1) = x(1);
Sy(N) = x(N);

i=2:N-1;
Sx(i) = -x(i);
Sy(i) = x(i);

dfdp = [ zeros(2*N,1) [Sx;Sy] ];
%%dfdp

% -----------------------------------------------------------------------

%%function df = fj(y,N,A,B,i_choose)
function df = fj(y,N,L,A,B,Dx,Dy,i_choose)

n = N*N;
x = y(1:n);
y = y(n+1:2*n);

%%Dx = 1;

%%L = 1;

%A = 1;
%%Dy = 0.02;

h  = 1/(N+1);
L2 = L*L;

cx = (Dx/L2)/(h*h);
cy = (Dy/L2)/(h*h);

%%h2 = h*h;
%%cx = Dx/h2;
%%cy = Dy/h2;

e  = ones(N,1);
TN = spdiags([e,-2*e,e], -1:1, N, N);
Tn = kron(eye(N),TN) + kron(TN,eye(N));
%%Tn = kron(eye(N),TN)/L2 + kron(TN,eye(N));

if i_choose == 1
  t = B*x + A*A*y + (B/A)*x.*x + 2*A*x.*y + x.*x.*y;

  dxdt = cx*Tn*x + t - x;
  dydt = cy*Tn*y - t;
  
  df = [dxdt; dydt];
  
  %fprintf('norm of df = %f\n',norm(df));

%%  norm(dxdt)
%%  norm(dydt)
%%pause
  
elseif i_choose == 2
  t_x = B + 2*(B/A)*x + 2*A*y + 2*x.*y;
  %t_x = -B + 2*x.*y; 
  t_y = A*A + 2*A*x + x.*x;
  %t_y = x.*x;
  
  f_reac_x = spdiags( t_x - 1, 0, n, n);
  f_reac_y = spdiags( t_y    , 0, n, n);
  g_reac_x = spdiags(-t_x    , 0, n, n);
  g_reac_y = spdiags(-t_y    , 0, n, n);

  dfdxy_reac = [f_reac_x f_reac_y
                g_reac_x g_reac_y];

  dfdxy_diff = blkdiag(cx*Tn, cy*Tn);
  df = dfdxy_diff + dfdxy_reac;
end
% --------------------------------------------------------------------------
function uf1 = user1(t,y,N,L,A,B,Dx,Dy)

uf1 = B-20.8; % step 1  BP1

function uf2 = user2(t,y,N,L,A,B,Dx,Dy)

uf2 = A-0.086; % step 2  BP2

function uf3 = user3(t,y,N,L,A,B,Dx,Dy)

uf3 = B-21.2;  % step 3  BP3

function uf4 = user4(t,y,N,L,A,B,Dx,Dy)

uf4 = A-0.089;  % step 4  BP4

% -------------------------------------------------------------------------
function normuser = usernorm(arg)

n = length(arg);
normuser = sqrt(2/n)*norm(arg);
