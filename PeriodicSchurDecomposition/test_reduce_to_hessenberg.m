clc;
m=9;
N=200;
G = zeros([N N m]);
for i=1:m
  G(:,:,i) = rand(N);
end
[H,Q]=reduce_to_hessenberg(G);


tolerance = eps(sqrt(N)*N);
maximum_error = 0;


for i=2:m
  assert(all(abs(Q(:,:,i)' * G(:,:,i) * Q(:,:,i-1)-H(:,:,i)) < tolerance,'all'))
  maximum_error_for_H_i = ...
    max(max(abs(Q(:,:,i)' * G(:,:,i) * Q(:,:,i-1)-H(:,:,i))));
  maximum_error = max(maximum_error_for_H_i, maximum_error); 
end
% note that Q(:,:,0) is defined as Q(:,:,m)
assert(all(abs(Q(:,:,1)' * G(:,:,1) * Q(:,:,m)-H(:,:,1)) < tolerance,'all'))
maximum_error_for_H_1 = max(max(abs(Q(:,:,1)' * G(:,:,1) * Q(:,:,m)-H(:,:,1))));
maximum_error = max(maximum_error_for_H_i, maximum_error); 
fprintf(['The maximum error in the periodic Hessenberg decomposition' ...
  ' in this test run is %.5e.\n'], maximum_error);
  