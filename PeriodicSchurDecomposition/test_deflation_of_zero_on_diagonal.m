clc;
m=5; % number of matrices
N=5; % size of matrix
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
% assert(all(abs(Q(:,:,1)' * H(:,:,1) * Q(:,:,m)-M(:,:,1)) < tolerance,'all'))
maximum_error_for_H_1 = max(max(abs(Q(:,:,1)' * H(:,:,1) * Q(:,:,m)-M(:,:,1))));
maximum_error = max(maximum_error_for_H_1, maximum_error); 
fprintf(['The maximum error in the deflation of a zero on the diagonal' ...
  ' in this test run is %.5e.\n'], maximum_error);

for i=1:m
  %max(max(Q(:,:,i)*Q(:,:,i)'-eye(size(Q,1))))
  assert(is_orthogonal(Q(:,:,i)))
end

% set entries that are very small to zero
% we hope that this way all entries that should be zero become zero
M(abs(M)<eps(N*N))=0;

for i=1:m-1
  % check if M_1 up to M_{m-1} are upper triangular
  assert(istriu(M(:,:,i)))
  % check if M_m is Hessenberg
  assert(bandwidth(M(:,:,m),'lower') == 1)
  % check if zeros have are on the correct place
  % on the subdiagonal of M_m
  assert(M(i_star  , i_star  , k_star) == 0)
  assert(M(i_star  , i_star-1, m     ) == 0)
  assert(M(i_star+1, i_star  , m     ) == 0)
end


function orthogonal=is_orthogonal(Q)
  N = size(Q,1);
  orthogonal = max(max(Q*Q'-eye(N))) < eps(sqrt(N)*N);
end

