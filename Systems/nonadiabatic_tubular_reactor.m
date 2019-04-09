function out = nonadiabatic_tubular_reactor
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
  out{10} =[]; %@usernorm
  out{11} = []; %@user function 1
end
% ----------------------------------------------------------------------
function du = fun_eval(~,x,D,P_em,P_eh,BETA,phi_0,GAMMA,B) % ignored parameter is t (time)
  N = length(x) / 2;
  h = 1 / (N-1);
  u = reshape(x,2,N);
  du = zeros(2,N);
  for i=1:N
    du(:,i) = reaction_term(u(:,i));
  end
  P = [P_em; P_eh];
  % at the left boundary, the term due to the derivative u_x is given by left
  % the Neumann BC:
  du(:,1) = du(:,1) + P .* (u(:,1) - 1);
  % at the right boundary, the term due to the derivative u_x is given by the
  % right Neumann BC, and equals zero.
  
  
  
  % At the left boundary the term due to diffusion term (the term with u_xx) is
  % computed using a virtual point u(:,0) such that u(:,2) == u(:,0) + 2 * h *
  % u_x(:,1). Hence u(:,1)_xx is approximated by:
  %
  %   (    u(:,0)                    - 2 * u(:,1) + u(:,2)) / h^2 
  % = (    u(:,2) - 2 * h * u_x(:,1) - 2 * u(:,1) + u(:,2)) / h^2
  % = (2 * u(:,2) - 2 * h * u_x(:,1) - 2 * u(:,1)         ) / h^2
  %
  % The Neumann BC u_x(:,1) is given as P .* (u(:,1) - 1). Hence u(:,1)_xx
  % equals:
  %
  % = (2 * u(:,2) - 2 * h * P .* (u(:,1) - 1) - 2 * u(:,1)) / h^2
  %
  du(:,1) = ...
    du(:,1) + (2 * u(:,2) - 2 * h * P .* (u(:,1) - 1) - 2 * u(:,1)) / h^2 ./ P;
  % At the right boundary the term due to u_xx is also computed using a virtual
  % point. The virtual point in this case is
  % u(:,N+1) == u(:,N-1) + 2 * h * u_x(:,N). But since the Neumann BC on the
  % right side is zero, the term with the Neumann BC vanishes:
  du(:,N) = ...
    du(:,N) + (2 * u(:,N-1)                           - 2 * u(:,N)) / h^2 ./ P;

  for i=2:N-1
    du(:,i) = du(:,i) + ...
         (u(:,i-1) - 2 * u(:,i) + u(:,i+1)) / (h^2) ./ P ... % diffusion term
       - (u(:,i+1)              - u(:,i-1)) /  h / 2;  % convection term
      
  end

  function f = reaction_term(u)
    y   = u(1);
    phi = u(2);
    D_y_exp_gamma_1_minus_1_over_phi = D * y * exp(GAMMA * (1 - 1/phi));
    f = [
                                 - D_y_exp_gamma_1_minus_1_over_phi;
      - BETA * (phi - phi_0) + B * D_y_exp_gamma_1_minus_1_over_phi
    ];
  end
  % transform du into column vector:
  du = du(:);
end
% --------------------------------------------------------------------------
function init(); end

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

