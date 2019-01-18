clc;
m=10;
N=10;
tolerance = eps(sqrt(N)*N);
G = zeros(N,N,m);
for i=1:m
  G(:,:,i) = rand(N);
end

[A,B,C] = mexSchurM(G);


