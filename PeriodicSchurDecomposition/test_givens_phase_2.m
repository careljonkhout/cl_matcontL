N = 4;

A = rand(N);

i = 4;
[c,s] = givens(A(i-1,i-1),A(i,i-1));
AA=givens_left(A,i-1,i,c,s)


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

function A = givens_left(A,i,k,c,s)
  for j=1:length(A)
    tau_1 = A(i,j);
    tau_2 = A(k,j);
    A(i,j) = c*tau_1 - s*tau_2;
    A(k,j) = s*tau_1 + c*tau_2;
  end
end
