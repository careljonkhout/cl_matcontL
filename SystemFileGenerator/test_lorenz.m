name = 'lorenz';
vars = 'x y z';
pars = 'sigma r b';
time = 't';
max_ord = 0;
rhs = {
  'sigma * (-x + y)'
  'r*x - y - x*z'
  '-b*z + x*y'
};

lorenz = System_of_ODEs.new(name,vars,pars,time,max_ord,rhs);
lorenz.generate_file()