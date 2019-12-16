function out = snode2d
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobian_params;
out{5} = @hessians;
out{6} = @hessians_params;
out{7} = @third_order_derivatives;
out{8} = @fourth_order_derivatives;
out{9} = @fifth_order_derivatives;
out{10} = [];% usernorm DV: this is the only difference between
% the standard matcont problem file and the cl_matcontL problem file
out{11}= @x1;
out{12}= @x2;
out{13}= @x3;
out{14}= @x4;
out{15}= @x5;

% Suppress function might not be used warning in this file:
%#ok<*DEFNU>
% Suppress Input argument might be unused warnings in this file: 
%#ok<*INUSD>
%#ok<*INUSL>

% --------------------------------------------------------------------------
function dydt = fun_eval(t,y, par_mu)
dydt =[
  par_mu - y(1)^2;
  - sin(y(2))];
% --------------------------------------------------------------------------
function options = init()
handles = feval(snode2d);
options = odeset(...
  'Jacobian', handles(3), 'JacobianP', handles(4), ...
  'Hessians', handles(5), 'HessiansP', handles(6));

% --------------------------------------------------------------------------
function jac = jacobian(t, y, par_mu)
jac = reshape([y(1).*-2.0,0.0,0.0,-cos(y(2))],[2,2]);
% --------------------------------------------------------------------------
function jacp=jacobian_params(t, y, par_mu)
jacp = [1.0;0.0];
% --------------------------------------------------------------------------
function hess = hessians(t, y, par_mu)
hess = reshape([-2.0,0.0,0.0,0.0,0.0,0.0,0.0,sin(y(2))],[2,2,2]);
% --------------------------------------------------------------------------
function hessp = hessians_params(t, y, par_mu)
hessp = [0.0;0.0];
%---------------------------------------------------------------------------
function d = third_order_derivatives(t, y, par_mu)
d = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,cos(y(2))],[2,2,2,2]);
%---------------------------------------------------------------------------
function d = fourth_order_derivatives(t, y, par_mu)
d = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-sin(y(2))],[2,2,2,2,2]);
%---------------------------------------------------------------------------
function d = fifth_order_derivatives(t, y, par_mu)
d = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-cos(y(2))],[2,2,2,2,2,2]);

function userfun1=x1(t, y, par_mu)
	userfun1=0;
function userfun2=x2(t, y, par_mu)
	userfun2=0;
function userfun3=x3(t, y, par_mu)
	userfun3=0;
function userfun4=x4(t, y, par_mu)
	userfun4=0;
function userfun5=x5(t, y, par_mu)
	userfun5=0;
