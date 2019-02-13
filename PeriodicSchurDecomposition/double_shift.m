function [M,Q] = double_shift(M,Q,L,U,a,b,c,d)
  N = size(M,1);
  m = size(M,3);
  n = U-L+1;
  H = eye(U-L+1);
  for i=1:m
    H = M(L:U,L:U,i)*H;
  end
  lambda_1 = a+d; %H(n-1,n-1) + H(n,n);
  lambda_2 = a*d - b*c; %H(n-1,n-1)*H(n,n) - H(n-1,n)*H(n,n-1);
  %H_bar = zeros(3,1);
  %H_bar(1) = ((H(n,n)-H(1,1))*(H(n-1,n-1)-H(1,1))-H(n,n-1)*H(n-1,n)) ... 
  %            / H(2,1)+H(1,2);
  %H_bar(2) =  H(2,2)-H(1,1)-(H(n,n)-H(1,1))-(H(n-1,n-1)-H(1,1));
  %H_bar(3) = H(2,1)*H(3,2);
  H_bar = (H-lambda_1 * eye(n))*(H-lambda_2 * eye(n));
  v = zeros(N,1);
  v(L:L+2) = H_bar(1:3,1);
  v(1) = v(1) + sign(v(1))*norm(v);
  v = v / norm(v)*sqrt(2);
  % Compute and apply shift
  M(:,:,m) = householder_left (v,M(:,:,m));
  M(:,:,1) = householder_right(v,M(:,:,1));
  Q(:,:,m) = householder_right(v,Q(:,:,m));
  % Restore upper Hessenberg structure
  for j=L:U-2
    for k=1:m-1
      % construct a Householder reflection to introduce zeros on positions 
      % (j+1,j) and (j+2,j) of M_k
      v = zeros(N,1);
      v(j:j+2) = M(j:j+2,j,k);
      v(j) = v(j) + sign(v(j))*norm(v);
      v = v / norm(v)*sqrt(2);    
      M(:,:,k  ) = householder_left (v,M(:,:,k  ));
      M(:,:,k+1) = householder_right(v,M(:,:,k+1));
      Q(:,:,k  ) = householder_right(v,Q(:,:,k  ));
    end
    if j < U-2
      % Construct a Householder transformation to introduce zeros at positions
      % (j+2,j) and (j+3,j) of M_m
      v = zeros(N,1);
      v(j+1:j+3) = M(j+1:j+3,j,m);
      v(j+1) = v(j+1) + sign(v(j+1))*norm(v);
      v = v / norm(v)*sqrt(2);    
      M(:,:,m) = householder_left (v,M(:,:,m));
      M(:,:,1) = householder_right(v,M(:,:,1));
      Q(:,:,m) = householder_right(v,Q(:,:,m));     
    else
      % Construct a Householder transformation to introduce a zero at position
      % (j+2,j) of M_m
      v = zeros(N,1);
      v(j+1:j+2) = M(j+1:j+2,j,m);
      v(j+1) = v(j+1) + sign(v(j+1))*norm(v);
      v = v / norm(v)*sqrt(2);    
      M(:,:,m) = householder_left (v,M(:,:,m));
      M(:,:,1) = householder_right(v,M(:,:,1));
      Q(:,:,m) = householder_right(v,Q(:,:,m));     
    end
  end


  %clc
  %M
  % Restore structure of column U-1

  for k=1:m-1
    % construct a Givens rotation to introduce a zero at position U,U-1 of M_k
    [c,s] = construct_givens(-M(U-1,U-1,k),M(U,U-1,k));
    M(:,:,k  ) = givens_left (M(:,:,k  ),U-1,U,c,s);
    M(:,:,k+1) = givens_right(M(:,:,k+1),U-1,U,c,s);
    Q(:,:,k  ) = givens_right(Q(:,:,k  ),U-1,U,c,s);
  end
  
  M = check_and_enforce_lower_triangular_and_hessenberg_structure(M);
end

function result = householder_left(v,M)
  result = M - v*(v'*M);
end

function result = householder_right(v,M)
  result = M - (M*v)*v';
end
