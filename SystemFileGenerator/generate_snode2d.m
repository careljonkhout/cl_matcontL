tic

name = 'snode2d';
vars = 'x y';
pars = 'mu';
time = 't';
max_ord = 5;

equations = {
    'x'' = mu - x^2'
    'y'' = - sin(y)'
};

 
snode2d_system = SystemFileGenerator.new(name, pars, time, max_ord, equations);
snode2d_system.generate_file
toc
