function out = Carel
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
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
function dydt = fun_eval(t,xxxxx, par_alpha, par_beta, par_mu, par_delta, par_gamma)
dydt=[par_mu - par_mu * xxxxx(1) - par_beta*(1 + par_delta * xxxxx(4)) * xxxxx(1) * xxxxx(3); par_beta*(1+par_delta*xxxxx(3)) *xxxxx(1)*xxxxx(3) - (par_mu + par_alpha)*xxxxx(2); par_alpha*xxxxx(2) - (par_mu+par_gamma)*xxxxx(3); xxxxx(4)-2*pi*xxxxx(5) - (xxxxx(4)^2+xxxxx(5)^2)*xxxxx(4); 2* pi * xxxxx(4) + xxxxx(5) - (xxxxx(4)^2+xxxxx(5)^2)*xxxxx(5)];
% --------------------------------------------------------------------------
function options = init()
handles = feval(Carel);
options = odeset(...
  'Jacobian', handles(3), 'JacobianP', handles(4), ...
  'Hessians', handles(5), 'HessiansP', handles(6));

% --------------------------------------------------------------------------
function jac = jacobian(t, xxxxx, par_alpha, par_beta, par_mu, par_delta, par_gamma)
jac = [];
% --------------------------------------------------------------------------
function jacp=jacobian_params(t, xxxxx, par_alpha, par_beta, par_mu, par_delta, par_gamma)
jacp = [];
% --------------------------------------------------------------------------
function hess = hessians(t, xxxxx, par_alpha, par_beta, par_mu, par_delta, par_gamma)
hess = [];
% --------------------------------------------------------------------------
function hessp = hessians_params(t, xxxxx, par_alpha, par_beta, par_mu, par_delta, par_gamma)
hessp = [];
%---------------------------------------------------------------------------
function d = third_order_derivatives(t, xxxxx, par_alpha, par_beta, par_mu, par_delta, par_gamma)
d = [];
%---------------------------------------------------------------------------
function d = fourth_order_derivatives(t, xxxxx, par_alpha, par_beta, par_mu, par_delta, par_gamma)
d = [];
%---------------------------------------------------------------------------
function d = fifth_order_derivatives(t, xxxxx, par_alpha, par_beta, par_mu, par_delta, par_gamma)
d = [];

function userfun1=x1(t, xxxxx, par_alpha, par_beta, par_mu, par_delta, par_gamma)
	userfun1=0;
function userfun2=x2(t, xxxxx, par_alpha, par_beta, par_mu, par_delta, par_gamma)
	userfun2=0;
function userfun3=x3(t, xxxxx, par_alpha, par_beta, par_mu, par_delta, par_gamma)
	userfun3=0;
function userfun4=x4(t, xxxxx, par_alpha, par_beta, par_mu, par_delta, par_gamma)
	userfun4=0;
function userfun5=x5(t, xxxxx, par_alpha, par_beta, par_mu, par_delta, par_gamma)
	userfun5=0;
