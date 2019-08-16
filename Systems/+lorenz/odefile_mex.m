function out = odefile_mex
out{1} = @init;
out{2} = @lorenz.dydt_mex;
out{3} = @lorenz.jacobian_mex;
out{4} = @lorenz.jacobian_params_mex;
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
function options = init()
handles = feval(lorenz);
options = odeset(...
  'Jacobian', handles(3), 'JacobianP', handles(4), ...
  'Hessians', handles(5), 'HessiansP', handles(6));



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
