function out = Thyroid_2d1
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
function dydt = fun_eval(t,y, par_k1, par_k2, par_sc1)
i_v1n = par_sc1 * 3.4e-12 * 1 * y(1) / (2.75 + 1*y(1));
i_v2n = 1e-10 * 1.5e-11 * y(2) / (1e-9 + 1.5e-11 * y(2));
dydt =[
  0.1*(( -i_v1n + par_k1 * i_v2n) + 5.7925 * 2.41692e2)/1 -       140 * y(1);
  0.4*(( -i_v2n + par_k2 * i_v1n) +       1e-18)/1.5e-11 - 1.000260643477599e+00 * 1.1e-6 * y(2)];
% --------------------------------------------------------------------------
function options = init()
handles = feval(Thyroid_2d1);
options = odeset(...
  'Jacobian', handles(3), 'JacobianP', handles(4), ...
  'Hessians', handles(5), 'HessiansP', handles(6));

% --------------------------------------------------------------------------
function jac = jacobian(t, y, par_k1, par_k2, par_sc1)
jac = reshape([(par_sc1.*-4.208996133570293e+15)./(y(1).*1.23794003928538e+28+3.404335108034796e+28)+par_sc1.*y(1).*1.0./(y(1).*1.23794003928538e+28+3.404335108034796e+28).^2.*5.210484838944022e+43-1.4e+2,(par_k2.*par_sc1.*4.208996133570293e+15)./(y(1).*4.642275147320176e+16+1.276625665513048e+17)-par_k2.*par_sc1.*y(1).*1.0./(y(1).*4.642275147320176e+16+1.276625665513048e+17).^2.*1.953931814604008e+32,(par_k1.*1.595073594941899e+15)./(y(2).*1.595073594941899e+26+1.063382396627933e+28)-par_k1.*y(2).*1.0./(y(2).*1.595073594941899e+26+1.063382396627933e+28).^2.*2.544259773280873e+41,y(2).*1.0./(y(2).*5.981525981032121e+14+3.987683987354748e+16).^2.*9.540974149803275e+29-1.595073594941899e+15./(y(2).*5.981525981032121e+14+3.987683987354748e+16)-1.100286707825359e-6],[2,2]);
% --------------------------------------------------------------------------
function jacp=jacobian_params(t, y, par_k1, par_k2, par_sc1)
jacp = reshape([(y(2).*1.595073594941899e+15)./(y(2).*1.595073594941899e+26+1.063382396627933e+28),0.0,0.0,(par_sc1.*y(1).*4.208996133570293e+15)./(y(1).*4.642275147320176e+16+1.276625665513048e+17),(y(1).*-4.208996133570293e+15)./(y(1).*1.23794003928538e+28+3.404335108034796e+28),(par_k2.*y(1).*4.208996133570293e+15)./(y(1).*4.642275147320176e+16+1.276625665513048e+17)],[2,3]);
% --------------------------------------------------------------------------
function hess = hessians(t, y, par_k1, par_k2, par_sc1)
hess = [];
% --------------------------------------------------------------------------
function hessp = hessians_params(t, y, par_k1, par_k2, par_sc1)
hessp = [];
%---------------------------------------------------------------------------
function d = third_order_derivatives(t, y, par_k1, par_k2, par_sc1)
d = [];
%---------------------------------------------------------------------------
function d = fourth_order_derivatives(t, y, par_k1, par_k2, par_sc1)
d = [];
%---------------------------------------------------------------------------
function d = fifth_order_derivatives(t, y, par_k1, par_k2, par_sc1)
d = [];

function userfun1=x1(t, y, par_k1, par_k2, par_sc1)
	userfun1=0;
function userfun2=x2(t, y, par_k1, par_k2, par_sc1)
	userfun2=0;
function userfun3=x3(t, y, par_k1, par_k2, par_sc1)
	userfun3=0;
function userfun4=x4(t, y, par_k1, par_k2, par_sc1)
	userfun4=0;
function userfun5=x5(t, y, par_k1, par_k2, par_sc1)
	userfun5=0;
