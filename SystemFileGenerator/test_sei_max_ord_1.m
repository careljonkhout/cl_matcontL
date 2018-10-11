tic

rhs=["mu - mu * s - beta*(1 + delta * u) * s * i"
    "beta*(1+delta*i) *s*i - (mu + alpha)*e"
    "alpha*e - (mu+gamma)*i"
    "u-2*pi*v - (u^2+v^2)*v"
    "2* pi * u + v - (u^2+v^2)*u"];
 
System_of_ODEs("Carel","s e i u v","alpha beta mu delta gamma","t",1,rhs)

toc

% max_order = 3 takes
% Elapsed time is 44.492268 seconds.