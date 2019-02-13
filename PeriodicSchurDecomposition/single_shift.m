function [M,Q]=single_shift(M,Q,mu,L,U)
  N = size(M,1);
  m = size(M,3);
  H = 1;
  for i=1:m-1
    H = H*M(L,L,i);
  end
  % Compute and apply shift
  H = M(L:L+1,L,m) * H;
  [c,s] = construct_givens(-H(1)+mu,H(2));
  M(:,:,m) = givens_left (M(:,:,m),L,L+1,c,s);
  M(:,:,1) = givens_right(M(:,:,1),L,L+1,c,s);
  Q(:,:,m) = givens_right(Q(:,:,m),L,L+1,c,s);
  % Restore upper Hessenberg structure
  
  for k=1:m-1
    % construct a Givens rotation to introduce a zero at position L+1,L of M_k
    [c,s] = construct_givens(-M(L,L,k),M(L+1,L,k));
    M(:,:,k  ) = givens_left (M(:,:,k  ),L,L+1,c,s);
    assert(abs(M(L+1,L,k))<eps(N));
    M(L+1,L,k) = 0;
    M(:,:,k+1) = givens_right(M(:,:,k+1),L,L+1,c,s);
    Q(:,:,k  ) = givens_right(Q(:,:,k  ),L,L+1,c,s);

  end
  
  for j=L:U-2
    % Restore structure of column j
    if j < U-1    
      % construct Householder reflection to introduce a zero at position
      % j+2,j of M_m
      x = M(:,j,m);       % x is the j-th column of M_m
      v = zeros(N,1);
      v(j+1:N) = x(j+1:N);
      v(j+1) = v(j+1) + sign(v(j+1))*norm(v);
      v = v / norm(v)*sqrt(2);
      % apply Householder reflection
      M(:,:,m) = householder_left( v,M(:,:,m));
      assert(M(j+2,j,m) < eps(N));
      M(j+2,j,m) = 0;    
      M(:,:,1) = householder_right(v,M(:,:,1));
      Q(:,:,m) = householder_right(v,Q(:,:,m));
      
    end

    for k=1:m-1
      % construct a Givens rotation to introduce a zero at position j+1,j of M_k
      [c,s] = construct_givens(-M(j+1,j+1,k),M(j+2,j+1,k));
      M(:,:,k  ) = givens_left (M(:,:,k  ),j+1,j+2,c,s);

      M(j+1,j,k) = 0; 
      M(:,:,k+1) = givens_right(M(:,:,k+1),j+1,j+2,c,s);
      assert(M(j+2,j,m) < eps(N));
      if j ==95 && ~ all(all(abs(tril(M(:,:,m),-2))<eps(m*N*N)))
        1;
      end
      Q(:,:,k  ) = givens_right(Q(:,:,k  ),j+1,j+2,c,s);

    end
  end
 % if ~ all(all(abs(tril(M(:,:,m),-2))<eps(m*N*N)))
  %  1;
  %end
  M = check_and_enforce_lower_triangular_and_hessenberg_structure(M);
end

function result = householder_left(v,M)
  result = M - v*(v'*M);
end

function result = householder_right(v,M)
  result = M - (M*v)*v';
end
