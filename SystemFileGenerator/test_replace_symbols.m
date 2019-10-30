clc
cd /home/carel/Documents/cl_matcontL/SystemFileGenerator
mex('replace_symbols.c', '-g');
% replace_symbols( ...
%  'a.^4.*y.^4.*sin(a.*yy(1).*y),b.^2.*y.^2.*sin(b.*y(1).^2.*y)', ...
%   {'x';'y';'z';'a'}, {'p_x', 'p_y', 'p_z', 'par_a'})
% replace_symbols( ...
%  'a', ...
%   {'x';'y';'z';'a'}, {'p_x', 'p_y', 'p_z', 'par_a'})
replace_symbols( ...
  'sin(a*x*x*y*aa*a*y)', ...
   {'x';'y';'z';'a'}, {'p_x', 'p_y', 'p_z', 'par_a'})