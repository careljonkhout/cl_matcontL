N = 3;

handles = fusion_symbolic();

L         = 10;
cn        = 1.1;
cT        = 0.9;
eps       = 0.05;
D0        = 1.9;
D1        = -1.1;
ZS        = 0;
Gamma_inf = -0.8;
lambda_n  = 1.25;
lambda_T  = 1.5;
gamma     = 1.6666666667;
mu        = 0.05;
zeta      = 1.1;
D2        = 0;
epsilon   = 0.05;

y = sym('y', [1 3*(N-1)]);
syms a b q_inf

dydt = feval(handles{2}, 0, y, N, ...
                 Gamma_inf, q_inf, D0, D1, D2, a, b, zeta, mu, ...
                 epsilon, ZS, gamma, lambda_n, lambda_T, cn, cT);
               
disp('simplifying dydt_sym')

dydt = simplify(dydt);
vars = char(y);
vars = vars(10:end-3);

tic
max_ord = 1;

name = sprintf('fusion_N_%d_max_ord_%d', N, max_ord);

pars = 'a b q_inf';
time = 't';

eqtns = cell(length(dydt),1);

for i = 1:length(dydt)
  eqtns{i} = [char(y(i)) '''=' char(dydt(i))];
end
output = 'cvode';

system = SystemFileGenerator.new( ...
                name, pars, time, max_ord, eqtns, output);
output = 'odefile';

toc
