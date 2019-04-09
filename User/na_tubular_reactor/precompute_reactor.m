% test of nonadiabatic_tubular_reactor odefile
clc
clear global
N = 100;                     
odefile = @nonadiabatic_tubular_reactor_symbolic;
D     = sym('D')    ;
P_em  = sym('P_em') ;
P_eh  = sym('P_eh') ;
BETA  = sym('BETA') ;
phi_0 = sym('phi_0');
GAMMA = sym('GAMMA');
B     = sym('B'    );

ode_parameters = {D; P_em; P_eh; BETA; phi_0; GAMMA; B};


handles = feval(odefile);
dydt = handles{2};
dydt_0 = feval(dydt,0,sym('u',[2*N 1]),ode_parameters{:})


