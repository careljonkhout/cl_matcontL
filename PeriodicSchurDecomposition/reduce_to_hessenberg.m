% reduces the array of matrices G to periodic upper Hessenberg form
%
% input: an array of m NxN nonsingular matrices G 
% output: arrays H and Q of m NxN matrices, such that the matrices in H are
% upper Hessenberg, and Q such that H(i) = Q(i)' * G(i) * Q(i-1) with 
% Q(0) := Q(m)
%
% based on algorithm 5.1 in paragraph 5.3.1 on page 219 of Kurt Lust' Phd
% thesis "Numerical bifurcation analysis of periodic solutions of partial
% differential equations" completed in 1997
function [H,Q]= reduce_to_hessenberg(G)
  N = size(G,1);
  m = size(G,3);
  H = G;
  Q = zeros(size(G));
  for i=1:m
    Q(:,:,i) = eye(N);
  end
  for i=1:N-1
    for j = 1:m-1
      % construct Householder reflection
      x = H(:,i,j);
      v = zeros(N,1);
      v(i:N) = x(i:N);
      v(i) = v(i) + sign(v(i)) * norm(v);
      v = v / norm(v);
      Householder_reflection = eye(N) - 2*(v*v');
      % apply Householder reflection
      H(:,:,j)   = Householder_reflection * H(:,:,j);
      H(:,:,j+1) = H(:,:,j+1) * Householder_reflection';
      Q(:,:,j)   = Q(:,:,j)   * Householder_reflection';
    end
    if i < N-1
      % construct Householder reflection      
      x = H(:,i,m);
      v = zeros(N,1);
      v(i+1:N) = x(i+1:N);
      v(i+1) = v(i+1) + sign(v(i+1)) * norm(v);
      v = v / norm(v);
      Householder_reflection = eye(N) - 2*(v*v');
      % apply Householder reflection
      H(:,:,m) = Householder_reflection * H(:,:,m);
      H(:,:,1) = H(:,:,1) * Householder_reflection';
      Q(:,:,m) = Q(:,:,m) * Householder_reflection';
    end
  end
end
  