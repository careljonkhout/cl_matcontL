handles = fusion;
dydt = handles{2};

syms a b q_inf

Gamma_inf = -0.8;
%q_inf     = -0.72;
D0        = 1.9;
D1        = -1.1;
D2        = 0;
%a         = -1;
%b         = -0.3;
zeta1     = 1.1;
mu1       = 0.05;
epsilon   = 0.05;
ZS        = 0;
gamma1    = 5/3;
lambda_n  = 1.25;
lambda_T  = 1.5;
c_n       = 1.1;
c_T       = 0.9;

N=50;

ode_parameters ={ 
  N, Gamma_inf, q_inf, D0, D1, D2, a, b, zeta1, mu1, epsilon, ZS, gamma1, ...
  lambda_n, lambda_T, c_n, c_T};

xxxxx = sym('xxxxx', [1 3*(N-1)]);


