clc;
m=5;
N=5;
tolerance = eps(sqrt(N)*N);
G = zeros([N N m]);
for i=1:m
  G(:,:,i) = rand(N);
end


[H,~]=reduce_to_hessenberg(G);

k_star=round(m/2);
i_star= round(N/2);

H(i_star,i_star,k_star) = 0;

Q = zeros([N N m]);
for i=1:m
  Q(:,:,i) = eye(N);
end



H(abs(H)<eps(100))=0;


[M,Q]=deflation_of_zero_on_diagonal(H,Q,k_star,i_star);


maximum_error = 0;
for i=2:m
  %assert(all(abs(Q(:,:,i)' * H(:,:,i) * Q(:,:,i-1)-M(:,:,i)) < tolerance,'all'))
  i;
  maximum_error_for_H_i = ...
    max(max(abs(Q(:,:,i)' * H(:,:,i) * Q(:,:,i-1)-M(:,:,i))));
  maximum_error = max(maximum_error_for_H_i, maximum_error); 
end
% note that Q(:,:,0) is defined as Q(:,:,m)
%assert(all(abs(Q(:,:,1)' * H(:,:,1) * Q(:,:,m)-M(:,:,1)) < tolerance,'all'))
%maximum_error_for_H_1 = max(max(abs(Q(:,:,1)' * H(:,:,1) * Q(:,:,m)-M(:,:,1))));
%maximum_error = max(maximum_error_for_H_1, maximum_error); 
fprintf(['The maximum error in the deflation of a zero on the diagonal' ...
  ' in this test run is %.5e.\n'], maximum_error);

for i=1:m
  %max(max(Q(:,:,i)*Q(:,:,i)'-eye(size(Q,1))))
  assert(is_orthogonal(Q(:,:,i)))
end

M(abs(M)<eps(100))=0;

function orthogonal=is_orthogonal(Q)
  N = size(Q,1);
  orthogonal = max(max(Q*Q'-eye(N))) < eps(sqrt(N)*N);
end

