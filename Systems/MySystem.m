function out = MySystem
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobian_params;
out{5} = @hessians;
out{6} = @hessians_params;
out{7} = @third_ord_derivatives;
out{8} = @fourth_ord_derivatives;
out{9} = @fifth_ord_derivatives;
out{10} = [];% usernorm DV: this is the only difference between
% the standard matcont problem file and the cl_matcontL problem file
out{11}= @x1;
out{12}= @x2;
out{13}= @x3;
out{14}= @x4;
out{15}= @x5;

% Suppress function might not be used warning in this file:
%#ok<*DEFNU>
% Suppress Input argument might be unused warning in this file: 
%#ok<*INUSD>

% --------------------------------------------------------------------------
function dydt = fun_eval(t,xxxxx, par_a, par_b)
dydt=[par_a*xxxxx(2), par_b*xxxxx(1)*(1-xxxxx(1))];
% --------------------------------------------------------------------------
%function [tspan,y0,options] = init
function y0 = init(par_alpha,par_beta)
%handles = feval(adapt22);
%y0=[0,0,0];
%options = odeset('Jacobian',handles(3),'JacobianP',handles(4),
%        'Hessians',handles(5),'HessiansP',handles(6));
%tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t, xxxxx, par_a, par_b)
jac=[[0, par_a], [- par_b*xxxxx(1) - par_b*(xxxxx(1) - 1), 0]];
% --------------------------------------------------------------------------
function jacp=jacobian_params(t, xxxxx, par_a, par_b)
jacp=[[xxxxx(2), 0], [0, -xxxxx(1)*(xxxxx(1) - 1)]];
% --------------------------------------------------------------------------
function hess = hessians(t, xxxxx, par_a, par_b)
hess=[[[0, 0], [0, 0]],[[-2*par_b, 0], [0, 0]]];
% --------------------------------------------------------------------------
function hessp = hessians_params(t, xxxxx, par_a, par_b)
hessp=[[[0, 0], [0, 0]],[[0, 0], [0, 0]]];
%---------------------------------------------------------------------------
function d = third_ord_derivatives(t, xxxxx, par_a, par_b)
d = [[[[0,0],[0,0]],[[0,0],[0,0]]],[[[0,0],[0,0]],[[0,0],[0,0]]]],;
%---------------------------------------------------------------------------
function d = fourth_ord_derivatives(t, xxxxx, par_a, par_b)
d = [[[[[0,0],[0,0]],[[0,0],[0,0]]],[[[0,0],[0,0]],[[0,0],[0,0]]]],[[[[0,0],[0,0]],[[0,0],[0,0]]],[[[0,0],[0,0]],[[0,0],[0,0]]]]],;
%---------------------------------------------------------------------------
function d = fifth_ord_derivatives(t, xxxxx, par_a, par_b)
d = [[[[[[0,0],[0,0]],[[0,0],[0,0]]],[[[0,0],[0,0]],[[0,0],[0,0]]]],[[[[0,0],[0,0]],[[0,0],[0,0]]],[[[0,0],[0,0]],[[0,0],[0,0]]]]],[[[[[0,0],[0,0]],[[0,0],[0,0]]],[[[0,0],[0,0]],[[0,0],[0,0]]]],[[[[0,0],[0,0]],[[0,0],[0,0]]],[[[0,0],[0,0]],[[0,0],[0,0]]]]]],;

function userfun1=x1(t, xxxxx, par_a, par_b)
	userfun1=0;
function userfun2=x2(t, xxxxx, par_a, par_b)
	userfun2=0;
function userfun3=x3(t, xxxxx, par_a, par_b)
	userfun3=0;
function userfun4=x4(t, xxxxx, par_a, par_b)
	userfun4=0;
function userfun5=x5(t, xxxxx, par_a, par_b)
	userfun5=0;
