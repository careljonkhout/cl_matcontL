function jac = contjac(x, varargin)
% Calculates jacobian matrix of F(x), which is to be found in curve file.
% The function handle to the function that defines the jacobian matrix of F(x),
% is stored in the field curve_jacobian of the global struct cds

global cds
jac =  feval(cds.curve_jacobian, x, varargin{:});