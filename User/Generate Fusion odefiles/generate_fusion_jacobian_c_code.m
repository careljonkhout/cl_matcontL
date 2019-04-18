N               = 50;
odefile         = str2func(sprintf('fusion_precomputed_with_sage_N_%d',N));
syms              a b q_inf
xxxxx           = sym('xxxxx', [1 3*(N-1)]);
handles             = feval(odefile);
jacobian_function   = handles{3};
jacobian = feval(jacobian_function,0,xxxxx,a,b,q_inf);
jacobian_c_code     = ccode(jacobian);

for i=3*(N-1):-1:1
  find    = sprintf('xxxxx%d',i);
  replace = sprintf('xxxxx[%d]',i-1);
  jacobian_c_code = strrep(jacobian_c_code, find, replace);
end
disp(jacobian_c_code(1:100))
pattern     = 'jacobian\[(\d+)\]\[(\d+)\]';
replacement = 'jacobian[$1 + $2 * 147]';
jacobian_c_code = regexprep(jacobian_c_code, pattern, replacement);
disp(jacobian_c_code(1:100))


filename        = sprintf('fusion_jacobian_c_code_N_%d',N);
fid             = fopen(filename,'w');
fprintf(fid,jacobian_c_code);
fclose(fid);
