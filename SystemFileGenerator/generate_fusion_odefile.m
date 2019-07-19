N = 25;

handles = fusion_symbolic();

dydt = handles{2};

load('fusion_standard_parameters')

L         = 10;
cn        = 1.1;
cT        = 0.9;
eps       = 0.05;
D0        = 1.9;
ZS        = 0;
Gamma_inf = -0.8;
lambda_n  = 1.25;
lambda_T  = 1.5;
gamma     = 1.6666666667;
mu        = 0.05;
zeta      = 1.1;
D2        = 0;

y = sym('y', [1 3*(N-1)]);
syms a b q_inf

dydt_sym = feval(handles{2}, 0, y, N, ...
                 Gamma_inf, q_inf, D0, D1, D2, a, b, zeta, mu, ...
                 epsilon, ZS, gamma, lambda_n, lambda_T, cn, cT);
               
disp('simplifying dydt_sym')
tic
oldVal = sympref('FloatingPointOutput',true);
dydt_sym = simplify(dydt_sym);
toc
%my_jacobian    = jacobian(dydt_sym, xxxxx);
%size(char(my_jacobian))
%simplified_jac = simplify(my_jacobian, 'Seconds', 2);
%size(char(simplified_jac))
%my_hessian     = hessian(dydt_sym, xxxxx);
vars = char(y);
vars = vars(10:end-3);

tic
max_ord = 1;

name = sprintf('fusion_N_%d_max_ord_%d', N, max_ord);

pars = 'a b q_inf';
time = 't';


rhs = cell(length(dydt_sym),1);

for i=1:length(dydt_sym)
  rhs{i} = char(dydt_sym(i));
end

c_code = false;

fusion_system = System_of_ODEs.new(name, vars, pars, time, max_ord, rhs, c_code);
fusion_system.generate_file()
toc
sympref('FloatingPointOutput',oldVal);
