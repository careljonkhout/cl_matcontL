N               = 75;
odefile         = str2func(sprintf('fusion_precomputed_with_sage_N_%d',N));
syms              a b q_inf
xxxxx           = sym('xxxxx', [1 3*(N-1)]);
handles         = feval(odefile);
dydt_function   = handles{2};
dydt = feval(dydt_function,0,xxxxx,a,b,q_inf);
dydt_c_code     = ccode(dydt);
disp(dydt_c_code(10000:10050))
for i=3*(N-1):-1:1
  dydt_c_code = ...
    strrep(dydt_c_code,sprintf('xxxxx%d',i),sprintf('xxxxx[%d]',i-1));
  dydt_c_code = ...
    strrep(dydt_c_code,sprintf('dydt[%d][0]',i-1),sprintf('dydt[%d]',i-1));
end
disp(dydt_c_code(10000:10050))


filename        = sprintf('fusion_dydt_ccode_N_%d',N);
fid             = fopen(filename,'w');
fprintf(fid,dydt_c_code);
fclose(fid);
