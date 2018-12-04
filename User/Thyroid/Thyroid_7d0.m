function out = Thyroid_7d0


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
function dydt = fun_eval(t,y,N,v0,v1,v01,v2,v3,v4,v5,v6,v7,b1,a1,Ss)

% Safety step which breaks simulation if negative concentration occurs
%if y(1)<0 || y(2)<0 || y(3)<0 || y(4)<0 || y(5)<0 test=y(6);
 % y
  %pause
%end

Thyroid_7d0_constants
 
t3n = y3*y(3)*1/(1+k31*IBS);
t3r = GR*t3n/(t3n+DR);
t4th= GT*y4*y(4)/(y4*y(4)+DT);
T4  = y1*y(1)*(1+k41*TBG+k42*TBPA);
%
f2coef = t4th*(y4*y(4)/(y4*y(4)+k));
TRH = v6*(1+v7*y(6))*TRH0;
f4f5coef = GH*v4*TRH/((v4*TRH+DH)*(1+Ss*y5*y(5)/(y5*y(5)+Ds))*(1+Ls*t3r));
%
f1 = (at*v0*v01*GT*(y4*y(4)/(y4*y(4)+DT))-bt*T4)/y1; % FT4
f2 = (a31*(v1*v01*GD1*(y1*y(1)/(y1*y(1)+km1)) + v1*v01*GD2*(y1*y(1)/(y1*y(1)+km2))...
    + v2*GT3*y4*y(4)/(y4*y(4)+DT)+ v2*GD1*f2coef/(km1+f2coef)...
    + v2*GD2*f2coef/(km2+f2coef))/(1+k30*TBG)-b31*y2*y(2))/y2;% FT3
f3 = (a32*v3*GD2*(y1*y(1)/(y1*y(1)+km2))-b32*v5*y3*y(3))/y3; % T3c
f4 = (as*f4f5coef - bs*y4*y(4))/y4;    % TSH
f5 = (as2*f4f5coef - bs2*y5*y(5))/y5;  % TSHz
%
f6 = a1*y(6) + b1*y(7)- y(6)*(y(6)^2 + y(7)^2);
f7 =-b1*y(6)+ a1*y(7) - y(7)*(y(6)^2 + y(7)^2);

%a2 = 0.5; a3 = 0.2;
%f6 = y(6)*(1-y(6)) - a1*y(6)*y(7)/(a3 + y(6));
%f7 = a2*y(7)*(1-y(7)/y(6));

dydt = [f1; f2; f3; f4; f5; f6; f7];

%y(6)
%y(7)
%dydt
%pause
% --------------------------------------------------------------------------
function y0 = init(N,v0,v1,v01,v2,v3,v4,v5,v6,v7,b1,a1,Ss)

y0 = zeros(N,1);

% Initial values

%       FT4          FT3          T3c        TSH         TSHz
%y0 = [1.0051e+00; 9.7762e-01; 9.8103e-01;  1.0067e+00; 8.7956e-01];
Thyroid_7d0_constants
%
y0 = y0_scaled;
% --------------------------------------------------------------------------
function dfdx = jacobian(t,y,N,v0,v1,v01,v2,v3,v4,v5,v6,v7,b1,a1,Ss)

% -------------------------------------------------------------------------
function dfdp = jacobianp(t,y,N,v0,v1,v01,v2,v3,v4,v5,v6,v7,b1,a1,Ss)

% -------------------------------------------------------------------------
function normuser = usernorm(arg)
n = length(arg);
normuser = sqrt(1/n)*norm(arg);
% -------------------------------------------------------------------------
function uf1 = user1(t,y,N,v0,v1,v01,v2,v3,v4,v5,v6,v7,b1,a1,Ss)
uf1 = L - 6.5e-002;
% ---------------------------------------------------------------------
% ----