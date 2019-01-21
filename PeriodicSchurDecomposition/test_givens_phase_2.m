N = 4;

A = rand(N);

i = 1;
[c,s] = construct_givens(-A(i,i),A(i+1,i));
AA=givens_left(A,i,i+1,c,s)



