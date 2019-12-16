tic

max_ord = 1;
name = 'sei';
pars = 'alpha beta mu delta gamma';
time = 't';


equations = {
%    'mu - mu * s - beta*(1 + delta * u) * s * i'
%    'beta*(1+delta*i) *s*i - (mu + alpha)*e'
%    'alpha*e - (mu+gamma)*i'
    's'' = mu * exp(-s) - mu - beta * (1 + delta * u) * exp(i)'
    'e'' = beta * (1 + delta * u) * exp(s + i - e) - mu - alpha'
    'i'' = alpha * exp(e - i) - mu - gamma'
    'u'' = u - 2*pi*v - (u^2 + v^2)*u'
    'v'' = 2*pi*u + v - (u^2 + v^2)*v'
};

 
sei_system = SystemFileGenerator.new(name, pars, time, max_ord, equations, ...
               'cvode');
sei_system.generate_file
toc

% max order 4 takes 72.990323 seconds.