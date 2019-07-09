% Odefile of 1-d Brusselator model
% This file uses a non-constant approximate solution as a starting guess.
function out = heterogeneous_x0
  out{1} = @init;
  out{2} = @Brusselator_1d.dydt;
  out{3} = @Brusselator_1d.jacobian;
  out{4} = []; %@Brusselator_1d.jacobian_params;
  out{5} = [];
  out{6} = [];
  out{7} = [];
  out{8} = [];
  out{9} = [];
  out{10} = @usernorm;
  out{11} = @user1;
end
% ------------------------------------------------------------------------------
function x0 = init(N,~,A,B,~,~) % unused argument are L, Dx, and Dy
  x0 = zeros(2*N,1);
  i=1:N;
  x0(i)   = A + 2*sin(pi*i/(N+1));
  x0(N+i) = B/A - 0.5*sin(pi*i/(N+1));
end
% ------------------------------------------------------------------------------
function normuser = usernorm(arg)
  n = length(arg);
  normuser = sqrt(1/n)*norm(arg);
end
% ------------------------------------------------------------------------------
function uf1 = user1(~,~,~,L,~,~,~,~) % unused arguments are t,x,N,A,B,Dx,Dy
  uf1 = L - 6.5e-002;
end
% ------------------------------------------------------------------------------