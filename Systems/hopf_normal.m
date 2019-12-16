function out = hopf_normal
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobian_params;
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
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
function dydt = fun_eval(t,y, par_b, par_sigma)
dydt =[
  par_b*y(1) - 2*pi*y(2) - par_sigma * (y(1)^2 + y(2)^2)*y(1);
  2*pi*y(1) + par_b*y(2) - par_sigma * (y(1)^2 + y(2)^2)*y(2)];
% --------------------------------------------------------------------------
function options = init()
handles = feval(hopf_normal);
options = odeset(...
  'Jacobian', handles(3), 'JacobianP', handles(4), ...
  'Hessians', handles(5), 'HessiansP', handles(6));

% --------------------------------------------------------------------------
function jac = jacobian(t, y, par_b, par_sigma)
jac = reshape([par_b-par_sigma.*y(1).^2.*2.0-par_sigma.*(y(1).^2+y(2).^2),pi.*2.0-par_sigma.*y(1).*y(2).*2.0,pi.*-2.0-par_sigma.*y(1).*y(2).*2.0,par_b-par_sigma.*y(2).^2.*2.0-par_sigma.*(y(1).^2+y(2).^2)],[2,2]);
% --------------------------------------------------------------------------
function jacp=jacobian_params(t, y, par_b, par_sigma)
jacp = reshape([y(1),y(2),-y(1).*(y(1).^2+y(2).^2),-y(2).*(y(1).^2+y(2).^2)],[2,2]);
% --------------------------------------------------------------------------
function hess = hessians(t, y, par_b, par_sigma)
hess = [];
% --------------------------------------------------------------------------
function hessp = hessians_params(t, y, par_b, par_sigma)
hessp = [];
%---------------------------------------------------------------------------
function d = third_order_derivatives(t, y, par_b, par_sigma)
d = [];
%---------------------------------------------------------------------------
function d = fourth_order_derivatives(t, y, par_b, par_sigma)
d = [];
%---------------------------------------------------------------------------
function d = fifth_order_derivatives(t, y, par_b, par_sigma)
d = [];

function userfun1=x1(t, y, par_b, par_sigma)
	userfun1=0;
function userfun2=x2(t, y, par_b, par_sigma)
	userfun2=0;
function userfun3=x3(t, y, par_b, par_sigma)
	userfun3=0;
function userfun4=x4(t, y, par_b, par_sigma)
	userfun4=0;
function userfun5=x5(t, y, par_b, par_sigma)
	userfun5=0;
