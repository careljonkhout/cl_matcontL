clc;
clear all; %#ok<CLALL>
m=4; % number of matrices
N=4; % size of matrix
tolerance = eps(sqrt(N)*N);
G = zeros([N N m]);
for i=1:m
  G(:,:,i) = rand(N);
end


[H,~]=reduce_to_hessenberg(G);


Q = zeros([N N m]);
for i=1:m
  Q(:,:,i) = eye(N);
end


prod = eye(N);
for i=m:-1:1
  prod = H(:,:,i) * prod;
end


evals = eig(prod);

H(abs(H)<eps(100))=0;
clc
M=H;
[M,Q]=double_shift(M,Q,1,N);

maximum_error = 0;
for i=2:m
  assert(...
    all(abs(Q(:,:,i)' * H(:,:,i) * Q(:,:,i-1)-M(:,:,i)) < tolerance,'all'),...
    ['Assertion failed: Q_%d * H_%d * Q_%d - M_%d exceeds tolerance' ...
     'max(max(abs(Q_%d * H_%d * Q_%d - M_%d))) = %.5f'], ...
     i,i,i-1,i,i,i,i-1, ...
     max(max(abs(Q(:,:,i)' * H(:,:,i) * Q(:,:,i-1)-M(:,:,i)))));
  maximum_error_for_H_i = ...
    max(max(abs(Q(:,:,i)' * H(:,:,i) * Q(:,:,i-1)-M(:,:,i))));
  maximum_error = max(maximum_error_for_H_i, maximum_error); 
end
% note that Q(:,:,0) is defined as Q(:,:,m)
assert(all(abs(Q(:,:,1)' * H(:,:,1) * Q(:,:,m)-M(:,:,1)) < tolerance,'all'), ...
  'Assertion failed: Q_1 * H_1 * Q_m - M_1 exceeds tolerance');
maximum_error_for_H_1 = max(max(abs(Q(:,:,1)' * H(:,:,1) * Q(:,:,m)-M(:,:,1))));
maximum_error = max(maximum_error_for_H_1, maximum_error); 
fprintf(['The maximum error in the deflation of a zero on the diagonal' ...
  ' in this test run is %.5e.\n'], maximum_error);

for i=1:m
  assert(is_orthogonal(Q(:,:,i)), ...
    ['Assertion failed: Q_%d is not orthogonal.' ...
     ' max(max(Q_%d * Q_%d - I)) = %.5e'], ...
       i,i,i, max(max(Q(:,:,i)*Q(:,:,i)'-eye(N))));
end

% set entries that are very small to zero
% we hope that this way all entries that should be zero become zero
M(abs(M)<eps(N*N))=0;

for i=1:m-1
  % check if M_1 up to M_{m-1} are upper triangular
  assert(istriu(M(:,:,i)),'Assertion failed: M_%d is not upper triangular',i);
end
% check if M_m is Hessenberg
assert(bandwidth(M(:,:,m),'lower') == 1, ...
  'Assertion failed: M_m is not Hessenberg');

function orthogonal=is_orthogonal(Q)
  N = size(Q,1);
  orthogonal = max(max(Q*Q'-eye(N))) < eps(N*N);
end
