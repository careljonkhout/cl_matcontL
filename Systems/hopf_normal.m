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
function dydt = fun_eval(t,xxxxx, par_b, par_sigma)
dydt=[par_b*xxxxx(1) - 2*pi*xxxxx(2) - par_sigma * (xxxxx(1)^2 + xxxxx(2)^2)*xxxxx(1); 2*pi*xxxxx(1) + par_b*xxxxx(2) - par_sigma * (xxxxx(1)^2 + xxxxx(2)^2)*xxxxx(2)];
% --------------------------------------------------------------------------
function options = init()
handles = feval(hopf_normal);
options = odeset(...
  'Jacobian', handles(3), 'JacobianP', handles(4), ...
  'Hessians', handles(5), 'HessiansP', handles(6));

% --------------------------------------------------------------------------
function jac = jacobian(t, xxxxx, par_b, par_sigma)
jac = reshape([[par_b - 2*par_sigma*xxxxx(1)^2 - par_sigma*(xxxxx(1)^2 + xxxxx(2)^2), - 2*pi - 2*par_sigma*xxxxx(1)*xxxxx(2)], [2*pi - 2*par_sigma*xxxxx(1)*xxxxx(2), par_b - 2*par_sigma*xxxxx(2)^2 - par_sigma*(xxxxx(1)^2 + xxxxx(2)^2)]],2,2)';
% --------------------------------------------------------------------------
function jacp=jacobian_params(t, xxxxx, par_b, par_sigma)
jacp = reshape([[xxxxx(1), -xxxxx(1)*(xxxxx(1)^2 + xxxxx(2)^2)], [xxxxx(2), -xxxxx(2)*(xxxxx(1)^2 + xxxxx(2)^2)]],2,2);
% --------------------------------------------------------------------------
function hess = hessians(t, xxxxx, par_b, par_sigma)
hess = [];
% --------------------------------------------------------------------------
function hessp = hessians_params(t, xxxxx, par_b, par_sigma)
hessp = [];
%---------------------------------------------------------------------------
function d = third_ord_derivatives(t, xxxxx, par_b, par_sigma)
d = [];
%---------------------------------------------------------------------------
function d = fourth_ord_derivatives(t, xxxxx, par_b, par_sigma)
d = [];
%---------------------------------------------------------------------------
function d = fifth_ord_derivatives(t, xxxxx, par_b, par_sigma)
d = [];

function userfun1=x1(t, xxxxx, par_b, par_sigma)
	userfun1=0;
function userfun2=x2(t, xxxxx, par_b, par_sigma)
	userfun2=0;
function userfun3=x3(t, xxxxx, par_b, par_sigma)
	userfun3=0;
function userfun4=x4(t, xxxxx, par_b, par_sigma)
	userfun4=0;
function userfun5=x5(t, xxxxx, par_b, par_sigma)
	userfun5=0;
