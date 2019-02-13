function [M, Q, is_deflated] = deflate_zeros(M,Q,is_deflated)
  k = 1; 
  N = size(M,1);
  m = size(M,3);
  while k <= m-1
    i = 1;
    while i <= N
      if ~ M(i,i,k) && ~ is_deflated(k,i) 
        is_deflated(k,i) = 1;
        [M,Q] = deflate_zero(M,Q,k,i);
        % start again from to beginning to look if any 
        % new zeros have been introduced on the diagonals 
        k = 0;
        break;
      end
      i = i + 1;
    end
    k = k + 1;
  end
end

function [M,Q]=deflate_zero(M,Q,k_star,i_star)
  N = size(M,1);
  m = size(M,3);
  for j=N-1:-1:i_star % j counting backwards
    for k=m:-1:k_star+1 % k counting backwards
      [c,s] = construct_givens(M(j+1,j+1,k),M(j+1,j,k));
      M(:,:,k)   = givens_right(M(:,:,k)  ,j,j+1,c,s);
      M(:,:,k-1) = givens_left (M(:,:,k-1),j,j+1,c,s);
      Q(:,:,k-1) = givens_right(Q(:,:,k-1),j,j+1,c,s);
    end
  end

  for j=N-1:-1:i_star+1  % j counting backwards
    for k=k_star:-1:1 % k counting backwards
      [c,s] = construct_givens(M(j+1,j+1,k),M(j+1,j,k));
      M(:,:,k) = givens_right(M(:,:,k),j,j+1,c,s);
      if k>1
        M(:,:,k-1) = givens_left (M(:,:,k-1),j,j+1,c,s);
        Q(:,:,k-1) = givens_right(Q(:,:,k-1),j,j+1,c,s);
      else
        M(:,:,m)   = givens_left (M(:,:,m),j,j+1,c,s);
        Q(:,:,m)   = givens_right(Q(:,:,m),j,j+1,c,s);
      end
    end
  end
 
  for i=2:i_star
     [c,s] = construct_givens(-M(i-1,i-1,m),M(i,i-1,m));
     M(:,:,m) = givens_left (M(:,:,m),i-1,i,c,s);
     M(:,:,1) = givens_right(M(:,:,1),i-1,i,c,s);
     Q(:,:,m) = givens_right(Q(:,:,m),i-1,i,c,s);
    for k=1:k_star-1
      [c,s] = construct_givens(-M(i-1,i-1,k),M(i,i-1,k));
      M(:,:,k  ) = givens_left (M(:,:,k  ),i-1,i,c,s);
      M(:,:,k+1) = givens_right(M(:,:,k+1),i-1,i,c,s);
      Q(:,:,k  ) = givens_right(Q(:,:,k  ),i-1,i,c,s);
    end
  end

  for i=2:i_star-1
    for k=k_star:m-1
      [c,s] = construct_givens(-M(i-1,i-1,k),M(i,i-1,k)); 
      M(:,:,k  ) = givens_left (M(:,:,k  ),i-1,i,c,s);
      M(:,:,k+1) = givens_right(M(:,:,k+1),i-1,i,c,s);
      Q(:,:,k  ) = givens_right(Q(:,:,k  ),i-1,i,c,s);
    end
  end
  M=check_and_enfore_lower_triangular_and_hessenberg_structure(M);
end
