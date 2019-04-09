% test of nonadiabatic_tubular_reactor odefile
clc
clear global
N = 5;                     
odefile = @nonadiabatic_tubular_reactor;
D     = 0;
P_em  = 5;
P_eh  = 5;
BETA  = 2.35;
phi_0 = 1;
GAMMA = 25;
B     = 0.5;

ode_parameters = {D; P_em; P_eh; BETA; phi_0; GAMMA; B};


handles = feval(odefile);
dydt = handles{2};
dydt_0 = feval(dydt,0,ones(2*N,1),ode_parameters{:});

assert(all(all(dydt_0 == zeros(2*N,1))));