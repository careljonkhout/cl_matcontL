tic

N = 5000;
y = sym('y', [1 2*N]);
syms L A B Dx Dy

dydt_sym = Brusselator_1d.dydt_reordered(0, y, N, L, A, B, Dx, Dy);
               
vars = char(y);
vars = vars(10:end-3);


max_ord = 1;

name = sprintf('brusselator_1d_N_%d', N);

pars = 'L A B Dx Dy';
time = 't';


rhs = arrayfun(@(e) char(e), dydt_sym, 'UniformOutput', false);

output = 'cvode';

system = System_of_ODEs.new(name, vars, pars, time, max_ord, rhs,output);
system.generate_file()
toc
