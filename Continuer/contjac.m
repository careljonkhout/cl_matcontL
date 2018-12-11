function jac = contjac(x, varargin)
%Calculates jacobian matrix of F(x), which is to be found in curve file, which is global

global cds
jac =  feval(cds.curve_jacobian, x, varargin{:});