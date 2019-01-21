N=5;
A=rand(N);

i=1;
k=3;

c = cos(1);
s = sin(1);


givens_mat = eye(N);
givens_mat(i,i) = c;
givens_mat(k,k) = c;
givens_mat(i,k) = s;
givens_mat(k,i) = -s;


assert(max(max(abs(A*givens_mat' - givens_right(A,i,k,c,s))))<eps(N));