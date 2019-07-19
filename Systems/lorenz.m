function out = lorenz
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobian_params;
out{5} = @hessians;
out{6} = @hessians_params;
out{7} = @third_order_derivatives;
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
function dydt = fun_eval(t,y, par_sigma, par_r, par_b)
dydt=[par_sigma * (-y(1) + y(2)); par_r*y(1) - y(2) - y(1)*y(3); -par_b*y(3) + y(1)*y(2)];
% --------------------------------------------------------------------------
function options = init()
handles = feval(lorenz);
options = odeset(...
  'Jacobian', handles(3), 'JacobianP', handles(4), ...
  'Hessians', handles(5), 'HessiansP', handles(6));

% --------------------------------------------------------------------------
function jac = jacobian(t, y, par_sigma, par_r, par_b)
jac = reshape([[-par_sigma, par_sigma, 0], [par_r - y(3), -1, -y(1)], [y(2), y(1), -par_b]],3,3);
% --------------------------------------------------------------------------
function jacp=jacobian_params(t, y, par_sigma, par_r, par_b)
jacp = reshape([[y(2) - y(1), 0, 0], [0, y(1), 0], [0, 0, -y(3)]],3,3);
% --------------------------------------------------------------------------
function hess = hessians(t, y, par_sigma, par_r, par_b)
hess = reshape([[[0, 0, 0], [0, 0, 0], [0, 0, 0]],[[0, 0, -1], [0, 0, 0], [-1, 0, 0]],[[0, 1, 0], [1, 0, 0], [0, 0, 0]]],3,3,3);
% --------------------------------------------------------------------------
function hessp = hessians_params(t, y, par_sigma, par_r, par_b)
hessp = reshape([[[0, 0, 0], [0, 0, 0], [0, 0, 0]],[[0, 0, 0], [0, 0, 0], [0, 0, 0]],[[0, 0, 0], [0, 0, 0], [0, 0, 0]]],3,3,3);
%---------------------------------------------------------------------------
function d = third_order_derivatives(t, y, par_sigma, par_r, par_b)
d = reshape([[[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]]],[[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]]],[[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]]]],3,3,3,3);
%---------------------------------------------------------------------------
function d = fourth_order_derivatives(t, y, par_sigma, par_r, par_b)
d = [];
%---------------------------------------------------------------------------
function d = fifth_order_derivatives(t, y, par_sigma, par_r, par_b)
d = [];

function userfun1=x1(t, y, par_sigma, par_r, par_b)
	userfun1=0;
function userfun2=x2(t, y, par_sigma, par_r, par_b)
	userfun2=0;
function userfun3=x3(t, y, par_sigma, par_r, par_b)
	userfun3=0;
function userfun4=x4(t, y, par_sigma, par_r, par_b)
	userfun4=0;
function userfun5=x5(t, y, par_sigma, par_r, par_b)
	userfun5=0;
