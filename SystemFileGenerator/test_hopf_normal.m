tic

name = 'hopf_normal';
pars = 'b sigma';
time = 't';
max_ord = 1;

equations ={
    'u'' = b*u - 2*pi*v - sigma * (u^2 + v^2)*u'
    'v'' = 2*pi*u + b*v - sigma * (u^2 + v^2)*v'
};

 
sei_system = SystemFileGenerator.new(name, pars, time, max_ord, equations);
sei_system.generate_file
toc

% max order 4 takes 72.990323 seconds.