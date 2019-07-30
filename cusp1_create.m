tic

name = 'cusp1';
vars = 'x';
pars = 'h r';
time = 't';
max_ord = 3;

rhs={
    'h + r*x - x^3'
};

cusp1 = System_of_ODEs.new(name, vars, pars, time, max_ord, rhs);
cusp1.generate_file
toc