name = 'Carel2';
s = System_of_ODEs.new(name,'x y','a b c','t',5,{'sin(a*x*y)','sin(x*x*y*b)'});
s.generate_file

fullpath = mfilename('fullpath');
my_path = fullpath(1:end-length(mfilename));

cd(my_path)
cd('..\Systems')
handles = eval(name);
rhs             = handles{2};
jacobian        = handles{3};
jacobian_params = handles{4};
hessians        = handles{5};
hessians_params = handles{6};
d3              = handles{7};
d4              = handles{8};
d5              = handles{9};

rhs_evaluated             = feval(rhs            ,0,zeros(2),1,1);
jacobian_evaluated        = feval(jacobian       ,0,zeros(2),1,1);
jacobian_params_evaluated = feval(jacobian_params,0,zeros(2),1,1);
hessians_evaluated        = feval(hessians       ,0,zeros(2),1,1);
hessians_params_evaluated = feval(hessians_params,0,zeros(2),1,1);
d3_evaluated              = feval(d3,             0,zeros(2),1,1);
d4_evaluated              = feval(d4,             0,zeros(2),1,1);
d5_evaluated              = feval(d5,             0,zeros(2),1,1);

assert(all(size(rhs_evaluated)             == [2,1]))
assert(all(size(jacobian_evaluated)        == [2,2]))
assert(all(size(jacobian_params_evaluated) == [2,3]))
assert(all(size(hessians_evaluated)        == [2,2,2]))
assert(all(size(hessians_params_evaluated) == [2,3,3]))
assert(all(size(d3_evaluated)              == [2,2,2,2]))
assert(all(size(d4_evaluated)              == [2,2,2,2,2]))
assert(all(size(d5_evaluated)              == [2,2,2,2,2,2]))