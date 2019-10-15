function out = odefile_mex
out{1} = @init;
out{2} = @reactor_N_50_cvode_max_ord_1.dydt_mex;
out{3} = @reactor_N_50_cvode_max_ord_1.jacobian_mex;
out{4} = @reactor_N_50_cvode_max_ord_1.jacobian_params_mex;
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
handles = feval(reactor_N_50_cvode_max_ord_1);
options = odeset(...
  'Jacobian', handles(3), 'JacobianP', handles(4), ...
  'Hessians', handles(5), 'HessiansP', handles(6));



function userfun1=x1(t, y, par_D, par_P_em, par_P_eh, par_BETA, par_phi_0, par_GAMMA, par_B)
	userfun1=0;
function userfun2=x2(t, y, par_D, par_P_em, par_P_eh, par_BETA, par_phi_0, par_GAMMA, par_B)
	userfun2=0;
function userfun3=x3(t, y, par_D, par_P_em, par_P_eh, par_BETA, par_phi_0, par_GAMMA, par_B)
	userfun3=0;
function userfun4=x4(t, y, par_D, par_P_em, par_P_eh, par_BETA, par_phi_0, par_GAMMA, par_B)
	userfun4=0;
function userfun5=x5(t, y, par_D, par_P_em, par_P_eh, par_BETA, par_phi_0, par_GAMMA, par_B)
	userfun5=0;
