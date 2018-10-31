%System_of_ODEs.new("Carel","x y","a b","t",4,["sin(a*x*y)","sin(x*x*y*b)"])
s = System_of_ODEs.new("Carel","x y","a b","t",4,["sin(a*x*y)","sin(x*x*y*b)"])
s.generate_file
