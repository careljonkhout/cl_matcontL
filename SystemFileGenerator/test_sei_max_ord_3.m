
max_ord = 3;
name = sprintf('SEI_max_ord_%d', max_ord);
parameters = "alpha beta mu delta gamma";
equations = [
    "s' = mu - mu * s - beta*(1 + delta * u) * s * i"
    "e' = beta*(1+delta*i) *s*i - (mu + alpha)*e"
    "i' = alpha*e - (mu+gamma)*i"
    "u' = u-2*pi*v - (u^2+v^2)*u"
    "v' = 2* pi * u + v - (u^2+v^2)*v"];
 
s = SystemFileGenerator.new(name, parameters, "t", max_ord, equations);
s.generate_file

old_path = pwd;

fullpath = mfilename('fullpath');
my_path = fullpath(1:end-length(mfilename));

cd(my_path)
cd('../Systems')
handles = eval(name);
rhs = handles{2};
jacobian = handles{3};

rhs_evaluated      = feval(rhs     ,0,zeros(5),1,1,1,1,1);
jacobian_evaluated = feval(jacobian,0,zeros(5),1,1,1,1,1);

assert(all(size(rhs_evaluated)      == [5,1]))
assert(all(size(jacobian_evaluated) == [5,5]))

cd(old_path)