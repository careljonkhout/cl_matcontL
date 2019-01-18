function [M,Q]=deflation_of_zero_on_diagonal(M,Q,k_star,i_star)
  N = size(M,1);
  m = size(M,3);
  for j=N-1:-1:i_star % j counting backwards
    for k=m:-1:k_star+1 % k counting backwards
      [c,s] = givens(-M(j+1,j+1,k),M(j+1,j,k));
      M(:,:,k)   = givens_right(M(:,:,k)  ,j,j+1,c,s);
      M(:,:,k-1) = givens_left (M(:,:,k-1),j,j+1,c,s);
      Q(:,:,k-1) = givens_right(Q(:,:,k-1),j,j+1,c,s);
    end
  end
  
  for j=N-1:-1:i_star+1
    for k=k_star:-1:1
      [c,s] = givens(-M(j+1,j+1,k),M(j+1,j,k));
      M(:,:,k) = givens_right(M(:,:,k),j,j+1,c,s);
      if k>1
        M(:,:,k-1) = givens_left(M(:,:,k-1),j,j+1,c,s);
        Q(:,:,k-1) = givens_right(Q(:,:,k-1),j,j+1,c,s);
      else
        M(:,:,m)   = givens_left(M(:,:,m),j,j+1,c,s);
        Q(:,:,m)   = givens_right(Q(:,:,m),j,j+1,c,s);
      end
    end
  end
  
 
  for i=2:i_star
    % something goes wrong here
        [c,s] = givens_lust(M(i-1,i-1,m),M(i,i-1,m));
     M(:,:,m) = givens_left_lust(M(:,:,m),i-1,i,c,s);
     M(:,:,1) = givens_right_lust(M(:,:,1),i-1,i,c,s);
     Q(:,:,m) = givens_right_lust(Q(:,:,m),i-1,i,c,s);

    for k=1:k_star-1
      
      %if k > 0; ke = k; else; ke = m; end
      [c,s] = givens_lust(M(i-1,i-1,k),M(i,i-1,k));
      M(:,:,k  ) = givens_left_lust (M(:,:,k  ),i-1,i,c,s);
      M(:,:,k+1) = givens_right_lust(M(:,:,k+1),i-1,i,c,s);
      Q(:,:,k  ) = givens_right_lust(Q(:,:,k  ),i-1,i,c,s);
    end
  
  end
  for i=2:i_star-1
    for k=k_star:m-1
      [c,s] = givens_lust(M(i-1,i-1,k),M(i,i-1,k));
      M(:,:,k  ) = givens_left_lust(M(:,:,k  ),i-1,i,c,s);
      M(:,:,k+1) = givens_right_lust(M(:,:,k+1),i-1,i,c,s);
      Q(:,:,k  ) = givens_right_lust(Q(:,:,k  ),i-1,i,c,s);
    end
  end
  
end

% Golub & van Loan
function [c,s] = givens(a,b)
  if b==0
    c=1;
    s=0;
  else
    if abs(b) > abs(a)
      t = -a/b;
      s = 1 / sqrt(1+t^2);
      c = s*t;
    else
      t = -b/a;
      c = 1 / sqrt(1+t^2);
      s = c*t;
    end
  end
end

function [c,s]=givens_lust(a,b)
  r = hypot(a,b);
  c = a/r;
  s = b/r;
end

function A = givens_left_lust(A,i,k,c,s)
  A=givens_left(A,i,k,c,-s);
  return
  givens_mat = eye(size(A));
  givens_mat(i,i) = c;
  givens_mat(k,k) = c;
  givens_mat(i,k) = s;
  givens_mat(k,i) = -s;
  A=givens_mat*A;
end

function A = givens_right_lust(A,i,k,c,s)
  givens_mat = eye(size(A));
  givens_mat(i,i) = c;
  givens_mat(k,k) = c;
  givens_mat(i,k) = s;
  givens_mat(k,i) = -s;
  A=A*givens_mat';
end


%
% Compute the matrix G * A
% where G is a givens rotation, i.e.
%
function A = givens_left(A,i,k,c,s)
  for j=1:length(A)
    tau_1 = A(i,j);
    tau_2 = A(k,j);
    A(i,j) = c*tau_1 - s*tau_2;
    A(k,j) = s*tau_1 + c*tau_2;
  end
end

%
% Compute the matrix  A * G'
%
%

function A = givens_right(A,i,k,c,s)
  for j=1:length(A)
    tau_1 = A(j,i);
    tau_2 = A(j,k);
    A(j,i) = c*tau_1 - s*tau_2;
    A(j,k) = s*tau_1 + c*tau_2;
  end
end

function orthogonal=isOrthogonal(Q)
  orthogonal = max(max(Q*Q'-eye(size(Q)))) < eps(size(Q,1)^2);
end
