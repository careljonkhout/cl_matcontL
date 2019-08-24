% Odefile of a spatial discretization of the PDE:
%
% u_t = u_xx + u - u^2
%
%
% The spatial domain is a one dimensional interval of length 1. The
% discretization uses finite differences on a equidistant mesh.

function out = scalar_reaction_diffusion
  out{1} = @init;
  out{2} = @d_u__d_t;
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
function du = d_u__d_t(~, u, L)
  % ignored parameter is t (time)
  N = length(u);
  
  u0 = 0;
  u1 = 0;
  du = zeros(N,1);
  
  for i=1:N
    du(i) = reaction_term(u(i));
  end

  du(1) = du(1) + (u0     - 2 * u(1) + u(2)) * (N+1)^2 / L / L;
  du(N) = du(N) + (u(N-1) - 2 * u(N) + u1  ) * (N+1)^2 / L / L;
  
  i = 2:N-1;
  du(i) = du(i) + (u(i-1) - 2 * u(i) + u(i+1)) * (N+1)^2 / L / L;

  function f = reaction_term(u)
    f = u - u^2;
  end
  % transform du into column vector:
end
% --------------------------------------------------------------------------
function init(); end

% --------------------------------------------------------------------------
function dfdx = jacobian(~, u, L) % unused argument is t
 N = length(u);

 c  = (N+1) * (N+1) / L / L;

 % A will store all the diagonals of the jacobian that have nonzero's.
 A = zeros(N,3);
 
 
 A(1:N, 1) =  c;
 A(1:N, 2) = - 2 * c + 1 - 2 * u;
 A(1:N, 3) =  c;
 
 % We compose the sparse matrix by passing the matrix A to the Matlab built-in
 % function spdiags. A now contains all the diagonals that have nonzero
 % elements.
 dfdx = spdiags(A, -1:1 , N, N);
end
% --------------------------------------------------------------------------

