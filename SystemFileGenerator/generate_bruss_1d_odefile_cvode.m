tic

N = 500;
y = sym('y', [1 2*N]);
syms L A B Dx Dy

dydt_sym = Brusselator_1d.dydt(0, y, N, L, A, B, Dx, Dy);
               
%disp('simplifying dydt_sym')
%tic
%oldVal = sympref('FloatingPointOutput',true);
%dydt_sym = simplify(dydt_sym);
%toc
vars = char(y);
vars = vars(10:end-3);


max_ord = 1;

name = sprintf('brusselator_1d_N_%d', N);

pars = 'L A B Dx Dy';
time = 't';


rhs = arrayfun(@(e) char(e), dydt_sym, 'UniformOutput', false);
%rhs = cell(length(dydt_sym),1);


%for i=1:length(dydt_sym)
%  rhs{i} = char(dydt_sym(i));
%end

output = 'cvode';

fusion_system = System_of_ODEs.new(name, vars, pars, time, max_ord, rhs,output);
fusion_system.generate_file()
toc
