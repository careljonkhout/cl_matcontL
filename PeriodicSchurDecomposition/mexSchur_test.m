clc;
m=10;
N=10;
tolerance = eps(sqrt(N)*N);
G = zeros(N,N,m);
for i=1:m
  G(:,:,i) = rand(N);
end


[H,Q]=reduce_to_hessenberg(G);
Hc = cell(m,1);
Qc = cell(m,1);

for i=1:m
  Hc{i} = H(:,:,i);
  Qc{i} = Q(:,:,i);
end

[A,B,C] = mexSchur(H,Q);


