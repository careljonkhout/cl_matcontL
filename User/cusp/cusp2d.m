function out = cusp2d
%
% Odefile of 2-d cusp model

out{1} = @init;
out{2} = @fun_eval;
out{3} = [];%@jacobian;
out{4} = [];%@jacobianp;
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
out{10}= [];@usernorm;
out{11}= [];@user1;
% ----------------------------------------------------------------------
function dydt = fun_eval(t,y,N,h,r)

% Safety step which breaks simulation if negative concentration occurs
%if y(1)<0 || y(2)<0 || y(3)<0 || y(4)<0 || y(5)<0 test=y(6);
 % y
  %pause
%end

f1 = h + r*y(1) - y(1)^3;
f2 = - 2*y(2);
dydt = [f1; f2];

%y(6)
%y(7)
%dydt
%pause
% --------------------------------------------------------------------------
function y0 = init(N,h,r)
% h = 0; r = 1;
y0 = zeros(N,1);
y0 = [1;0]
% --------------------------------------------------------------------------
function dfdx = jacobian(t,y,N,h,r)

% -------------------------------------------------------------------------
function dfdp = jacobianp(t,y,N,h,r)

% -------------------------------------------------------------------------
function normuser = usernorm(arg)
n = length(arg);
normuser = sqrt(1/n)*norm(arg);
% -------------------------------------------------------------------------
function uf1 = user1(t,y,N,h,r)
uf1 = L - 6.5e-002;
% ---------------------------------------------------------------------
% ----