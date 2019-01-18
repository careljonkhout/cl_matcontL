N = 5;

A = rand(N);

i = 1;
k = 2;

[c,s] = givens(A(k,k),A(k,i));

A=givens_right(A,i,k,c,s)


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


function A = givens_left_lust(A,i,k,c,s)
  A=givens_left(A,i,k,c,-s);
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
    A(i,j) =   c*tau_1 + s*tau_2;
    A(k,j) = - s*tau_1 + c*tau_2;
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
    A(j,i) =   c*tau_1 + s*tau_2;
    A(j,k) = - s*tau_1 + c*tau_2;
  end
end

