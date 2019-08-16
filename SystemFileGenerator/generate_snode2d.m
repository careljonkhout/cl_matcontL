tic

name = 'snode2d';
vars = 'x y';
pars = 'mu';
time = 't';
max_ord = 5;

rhs={
    'mu - x^2'
    '- sin(y)'
};

 
snode2d_system = System_of_ODEs.new(name, vars, pars, time, max_ord, rhs);
snode2d_system.generate_file
toc
