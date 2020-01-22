function out = respiration_161_mex
out{1} = @init;
out{2} = @respiration_161_dydt;
out{3} = @respiration_161_jacobian;
out{4} = @respiration_161_jacobian_params;
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
handles = feval(brusselator_1d_N_5);
options = odeset(...
  'Jacobian', handles(3), 'JacobianP', handles(4), ...
  'Hessians', handles(5), 'HessiansP', handles(6));



function userfun1=x1(t, y, par_L, par_A, par_B, par_Dx, par_Dy)
	userfun1=0;
function userfun2=x2(t, y, par_L, par_A, par_B, par_Dx, par_Dy)
	userfun2=0;
function userfun3=x3(t, y, par_L, par_A, par_B, par_Dx, par_Dy)
	userfun3=0;
function userfun4=x4(t, y, par_L, par_A, par_B, par_Dx, par_Dy)
	userfun4=0;
function userfun5=x5(t, y, par_L, par_A, par_B, par_Dx, par_Dy)
	userfun5=0;
