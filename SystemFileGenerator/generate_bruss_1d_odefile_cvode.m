tic

N = 10;
y = sym('y', [1 2*N]);
syms L A B Dx Dy

dydt = Brusselator_1d.dydt_reordered(0, y, N, L, A, B, Dx, Dy);
               
vars = char(y);
vars = vars(10:end-3);


max_ord = 1;

name = sprintf('brusselator_1d_N_%d', N);

pars = 'L A B Dx Dy';
time = 't';

eqtns = cell(length(dydt),1);

for i = 1:length(dydt)
  eqtns{i} = [char(y(i)) '''=' char(dydt(i))];
end
output = 'cvode';

system = SystemFileGenerator.new( ...
                name, pars, time, max_ord, eqtns, output);
system.generate_file()
toc
