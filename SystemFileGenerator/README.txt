tested only on matlab 2017a
  
To generate a matcont system file, start the GUI with the command

SystemsGUI.new
  
or if you want the use the command line interface execute the commands:

s = System_of_ODEs.new(<name>, <variables_str>, <parameters_str>, <time>, <maxOrder>, <rhs>)
s.generate_file
 
where <name> is the name of the system, which will be used for the filename. The matcont file generated represents the first order system of ODEs

dx/dt = f(t,x,a)

where f is a function from R^n to R^n and where n is the length of x

- f is specified by the variable rhs
- rhs should be supplied as an n by 1 string array, a 1 by n string array or an n by m char array, where m is the length of the longest element of rhs. 
- The elements of f should be valid matlab expressions, using only variables and parameters specified in variables_str and parameters_str.
- The n-th element of rhs represents the n-th element of f(x,a). THIS IS AN IMPORTANT DIFFERENCE FROM THE MATCONT SYSTEM EDITOR
- x is specified by the variable variables_str.
- a is specified by the variable parameters_str.
- variables_str and parameters_str should be supplied as a char array or a string, such that the variables/parameters are separated by spaces or comma's. Additional spaces or comma's between the 
variables/parameters are ignored.
- the symbol that represents t in dx/dt = f(t,x,a) is specified by the variable time.