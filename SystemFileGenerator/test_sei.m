tic

max_ord = 4;
name = sprintf('SEI_max_ord_%d', max_ord);
vars = 's e i u v';
pars = 'alpha beta mu delta gamma';
time = 't';


rhs={
%    'mu - mu * s - beta*(1 + delta * u) * s * i'
%    'beta*(1+delta*i) *s*i - (mu + alpha)*e'
%    'alpha*e - (mu+gamma)*i'
    'mu * exp(-s) - mu - beta * (1 + delta * u) * exp(i)'
    'beta * (1 + delta * u) * exp(s + i - e) - mu - alpha'
    'alpha * exp(e - i) - mu - gamma'
    'u - 2*pi*v - (u^2 + v^2)*u'
    '2*pi*u + v - (u^2 + v^2)*v'
};

 
sei_system = System_of_ODEs.new(name, vars, pars, time, max_ord, rhs);
sei_system.generate_file
toc

% max order 4 takes 72.990323 seconds.