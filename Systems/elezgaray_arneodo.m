function out = elezgaray_arneodo
%
% Odefile of a space discretization of a 1-d nonadiabatic tubular reactor.
% Discretized using a finite differences on a equidistant mesh.

  out{1} = @init;
  out{2} = @d_y__d_t;
  out{3} = @jacobian;
  out{4} = [];%@jacobianp;
  out{5} = [];
  out{6} = [];
  out{7} = [];
  out{8} = [];
  out{9} = [];
  out{10} =[]; %@usernorm
  out{11} = []; %@user function 1
end
% ----------------------------------------------------------------------
function dy = d_y__d_t(~, y, D, eps, alpha)
  % ignored parameter is t (time)
  N = length(y) / 2;
  h = 1 / (N+1);
  y0 = [-2 -4];
  y1 = [-2 -4];
  y = reshape(y,N,2);
  dy = zeros(N,2);
  
  for i=1:N
    dy(i,:) = reaction_term(y(i,:));
  end

  dy(1,:) = dy(1,:) + D * (y0       - 2 * y(1,:) + y(2,:))/ h / h;
  dy(N,:) = dy(N,:) + D * (y(N-1,:) - 2 * y(N,:) + y1    )/ h / h;
  
  
  i = 2:N-1;
  dy(i,:) = dy(i,:) + D * (y(i-1,:) - 2 * y(i,:) + y(i+1,:)) / (h^2);

  function f = reaction_term(y)
    u = y(1);
    v = y(2);
    
    f = [
      (v - (u^2+u^3)) / eps
      alpha - u
    ];
  end
  % transform du into column vector:
  dy = dy(:);
end
% --------------------------------------------------------------------------
function init(); end

% --------------------------------------------------------------------------
function dfdx = jacobian(~, y, D, eps, ~) % unused arguments are t and alpha
 N = length(y)/2;
 u = y(1:N)';
 
 h = 1 / (N+1) ;
 c  = D/(h*h);

 % Sparse jacobian
 A = zeros(2*N,3);
 
 
 % The jacobian is can be divided in four blocks of size N by N
 % the upper left  block corresponds to d(u_t) / d u
 % the upper right block corresponds to d(u_t) / d v
 % the lower left  block corresopnds to d(v_t) / d u
 % the lower right block corresponds to d(v_t) / d v
 
 % the upper left block contains terms due to the discretized diffusion term
 % u_xx, and a term to the reaction term:
 A(1:N-1,   2) =  c;
 A(1:N,     3) = - 2 * c - (2 * u  + 3 * u .* u ) / eps;
 A(1:N,     4) =  c;
 
 % the lower right block contains terms due to the discretized diffusion term:
 A(N+1:2*N, 2) =  c;
 A(N+1:2*N, 3) = -2*c;
 A(N+2:2*N, 4) =  c;
 
 % the lower left block contains a diagonal due to the reaction term:
 A(1:N,     1) = -1 * ones(N,1);
 % the upper right block contains a diagonal due to the reaction term
 A(N+1:2*N, 5) = 1 / eps;
 % We compose the sparse matrix by passing the matrix A to the Matlab built-in
 % function spdiags. A now contains all the diagonals that have nonzero
 % elements.
 dfdx = spdiags(A, [-N,-1:1,N] , 2*N, 2*N);
end
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
end
% -------------------------------------------------------------------------

