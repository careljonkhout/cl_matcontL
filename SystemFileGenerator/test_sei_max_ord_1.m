
max_ord = 1;
name = sprintf('sei_max_ord_%d', max_ord);

rhs=["mu - mu * s - beta*(1 + delta * u) * s * i"
    "beta*(1+delta*i) *s*i - (mu + alpha)*e"
    "alpha*e - (mu+gamma)*i"
    "u-2*pi*v - (u^2+v^2)*u"
    "2* pi * u + v - (u^2+v^2)*v"];
 
s = System_of_ODEs.new(name,"s e i u v","alpha beta mu delta gamma","t",1,rhs)
s.generate_file

old_path = pwd;

fullpath = mfilename('fullpath');
my_path = fullpath(1:end-length(mfilename));

cd(my_path)
cd('..\Systems')
handles = eval(name);
rhs = handles{2};
jacobian = handles{3};

rhs_evaluated      = feval(rhs     ,0,zeros(5),1,1,1,1,1);
jacobian_evaluated = feval(jacobian,0,zeros(5),1,1,1,1,1);

assert(all(size(rhs_evaluated)      == [5,1]))
assert(all(size(jacobian_evaluated) == [5,5]))

cd(old_path)
% max_order = 3 takes
% Elapsed time is 44.492268 seconds.