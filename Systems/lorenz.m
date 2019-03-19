function out = lorenz
out{1} = [];
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
function dydt = fun_eval(t,xxxxx, par_sigma, par_r, par_b)
dydt=[par_sigma * (-xxxxx(1) + xxxxx(2)); par_r*xxxxx(1) - xxxxx(2) - xxxxx(1)*xxxxx(3); -par_b*xxxxx(3) + xxxxx(1)*xxxxx(2)];
% --------------------------------------------------------------------------

% --------------------------------------------------------------------------
function jac = jacobian(t, xxxxx, par_sigma, par_r, par_b)
jac = [];
% --------------------------------------------------------------------------
function jacp=jacobian_params(t, xxxxx, par_sigma, par_r, par_b)
jacp = [];
% --------------------------------------------------------------------------
function hess = hessians(t, xxxxx, par_sigma, par_r, par_b)
hess = [];
% --------------------------------------------------------------------------
function hessp = hessians_params(t, xxxxx, par_sigma, par_r, par_b)
hessp = [];
%---------------------------------------------------------------------------
function d = third_ord_derivatives(t, xxxxx, par_sigma, par_r, par_b)
d = [];
%---------------------------------------------------------------------------
function d = fourth_ord_derivatives(t, xxxxx, par_sigma, par_r, par_b)
d = [];
%---------------------------------------------------------------------------
function d = fifth_ord_derivatives(t, xxxxx, par_sigma, par_r, par_b)
d = [];

function userfun1=x1(t, xxxxx, par_sigma, par_r, par_b)
	userfun1=0;
function userfun2=x2(t, xxxxx, par_sigma, par_r, par_b)
	userfun2=0;
function userfun3=x3(t, xxxxx, par_sigma, par_r, par_b)
	userfun3=0;
function userfun4=x4(t, xxxxx, par_sigma, par_r, par_b)
	userfun4=0;
function userfun5=x5(t, xxxxx, par_sigma, par_r, par_b)
	userfun5=0;
