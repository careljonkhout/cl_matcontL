N = 50;

handles = nonadiabatic_tubular_reactor_symbolic();

dydt = handles{2};

y = sym('y', [1 2*N]);
syms D P_em P_eh BETA phi_0 GAMMA_ B

dydt_sym = feval(handles{2}, 0, y, D,P_em,P_eh,BETA,phi_0,GAMMA_,B);
               
disp('simplifying dydt_sym')
tic
oldVal = sympref('FloatingPointOutput',true);
dydt_sym = simplify(dydt_sym);
toc
vars = char(y);
vars = vars(10:end-3);

tic
max_ord = 1;

name = sprintf('reactor_N_%d_cvode_max_ord_%d', N, max_ord);

pars = 'D P_em P_eh BETA phi_0 GAMMA_ B';
time = 't';


equations = cell(length(dydt_sym),1);

for i = 1:length(dydt_sym)
  equations{i} = [char(y(i)) '''=' char(dydt_sym(i))];
end

output = 'cvode';

system = SystemFileGenerator.new(name, pars, time, max_ord, equations, output);
system.generate_file()
toc
sympref('FloatingPointOutput',oldVal);
