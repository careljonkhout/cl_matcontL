N = 25;

odefile = str2func(['fusion_precomputed_with_sage_N_' num2str(N)]);

handles = feval(odefile);

dydt = handles{2};


xxxxx = sym('xxxxx', [1 3*(N-1)]);
syms a b q_inf

dydt_sym = feval(handles{2}, 0, xxxxx, a, b, q_inf);
disp('simplifying dydt_sym')
tic
dydt_sym = simplify(dydt_sym);
toc
%my_jacobian    = jacobian(dydt_sym, xxxxx);
%size(char(my_jacobian))
%simplified_jac = simplify(my_jacobian, 'Seconds', 2);
%size(char(simplified_jac))
%my_hessian     = hessian(dydt_sym, xxxxx);
vars = char(xxxxx);
vars = vars(10:end-3);

tic

name = 'fusion_N_25_with_2_ord_ders';

pars = 'a b q_inf';
time = 't';
max_ord = 2;

rhs = cell(length(dydt_sym),1);

for i=1:length(dydt_sym)
  rhs{i} = char(dydt_sym(i));
end



 
fusion_system = System_of_ODEs.new(name, vars, pars, time, max_ord, rhs);
fusion_system.generate_file()
toc