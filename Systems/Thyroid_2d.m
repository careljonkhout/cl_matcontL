function out = Thyroid_2d
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
function dydt = fun_eval(t,y, par_V2, par_s1, par_e1)
dydt=[(par_s1 - (205517389334486962890625*y(1)^2)/(604462909807314587353088*(y(1)^2 + 11/4)) + (100000000000*par_V2*y(2)^2)/(y(2)^2 + 1/1000000000))/y(1) - par_e1*y(1); ((4208996133570293*y(1)^2)/(3094850098213450687247810560*(y(1)^2 + 11/4)) - (2*par_V2*y(2)^2)/(5*(y(2)^2 + 1/1000000000)))/y(2)];
% --------------------------------------------------------------------------
function options = init()
handles = feval(Thyroid_2d);
options = odeset(...
  'Jacobian', handles(3), 'JacobianP', handles(4), ...
  'Hessians', handles(5), 'HessiansP', handles(6));

% --------------------------------------------------------------------------
function jac = jacobian(t, y, par_V2, par_s1, par_e1)
jac = reshape([-par_e1-1.0./y(1).^2.*(par_s1-(y(1).^2.*2.05517389334487e+23)./(y(1).^2.*6.044629098073146e+23+1.662273001970115e+24)+(par_V2.*y(2).^2.*1.0e+11)./(y(2).^2+1.0e-9))+(y(1).^3.*1.0./(y(1).^2.*6.044629098073146e+23+1.662273001970115e+24).^2.*2.484552783462535e+47-(y(1).*4.110347786689739e+23)./(y(1).^2.*6.044629098073146e+23+1.662273001970115e+24))./y(1),-(y(1).^3.*1.0./(y(1).^2.*3.094850098213451e+27+8.510837770086989e+27).^2.*2.605242419472011e+43-(y(1).*8.417992267140586e+15)./(y(1).^2.*3.094850098213451e+27+8.510837770086989e+27))./y(2),((par_V2.*y(2).*2.0e+11)./(y(2).^2+1.0e-9)-par_V2.*y(2).^3.*1.0./(y(2).^2+1.0e-9).^2.*2.0e+11)./y(1),-((par_V2.*y(2).*4.0)./(y(2).^2.*5.0+5.0e-9)-par_V2.*y(2).^3.*1.0./(y(2).^2.*5.0+5.0e-9).^2.*2.0e+1)./y(2)-1.0./y(2).^2.*((y(1).^2.*4.208996133570293e+15)./(y(1).^2.*3.094850098213451e+27+8.510837770086989e+27)-(par_V2.*y(2).^2.*2.0)./(y(2).^2.*5.0+5.0e-9))],[2,2]);
% --------------------------------------------------------------------------
function jacp=jacobian_params(t, y, par_V2, par_s1, par_e1)
jacp = reshape([(y(2).^2.*1.0e+11)./(y(1).*(y(2).^2+1.0e-9)),(y(2).*-2.0)./(y(2).^2.*5.0+5.0e-9),1.0./y(1),0.0,-y(1),0.0],[2,3]);
% --------------------------------------------------------------------------
function hess = hessians(t, y, par_V2, par_s1, par_e1)
hess = [];
% --------------------------------------------------------------------------
function hessp = hessians_params(t, y, par_V2, par_s1, par_e1)
hessp = [];
%---------------------------------------------------------------------------
function d = third_order_derivatives(t, y, par_V2, par_s1, par_e1)
d = [];
%---------------------------------------------------------------------------
function d = fourth_order_derivatives(t, y, par_V2, par_s1, par_e1)
d = [];
%---------------------------------------------------------------------------
function d = fifth_order_derivatives(t, y, par_V2, par_s1, par_e1)
d = [];

function userfun1=x1(t, y, par_V2, par_s1, par_e1)
	userfun1=0;
function userfun2=x2(t, y, par_V2, par_s1, par_e1)
	userfun2=0;
function userfun3=x3(t, y, par_V2, par_s1, par_e1)
	userfun3=0;
function userfun4=x4(t, y, par_V2, par_s1, par_e1)
	userfun4=0;
function userfun5=x5(t, y, par_V2, par_s1, par_e1)
	userfun5=0;
