tic

rhs=["mu - mu * s - beta*(1 + delta * u) * s * i"
    "beta*(1+delta*i) *s*i - (mu + alpha)*e"
    "alpha*e - (mu+gamma)*i"
    "u-2*pi*v - (u^2+v^2)*u"
    "2* pi * u + v - (u^2+v^2)*v"];
 
sei_system=System_of_ODEs.new("Carel","s e i u v", "alpha beta mu delta gamma","t",0,rhs)
sei_system.generate_file
toc

% max order 4 takes 72.990323 seconds.