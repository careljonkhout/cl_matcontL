tic

name = 'hopf_normal';
vars = 'u v';
pars = 'b sigma';
time = 't';
max_ord = 1;

rhs={
    'b*u - 2*pi*v - sigma * (u^2 + v^2)*u'
    '2*pi*u + b*v - sigma * (u^2 + v^2)*v'
};

 
sei_system = System_of_ODEs.new(name, vars, pars, time, max_ord, rhs);
sei_system.generate_file
toc

% max order 4 takes 72.990323 seconds.