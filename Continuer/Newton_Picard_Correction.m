function pout = Newton_Picard_Correction(x0)

global cds

% 1

phi = feval(cds.curve_func, x)

jac =  feval(cds.curve_jacobian, x, varargin{:});
