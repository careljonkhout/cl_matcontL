N = 5;

A = rand(N);

i = 1;
k = 2;

[c,s] = givens(-A(k,k),A(k,i));

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

function A = givens_right(A,i,k,c,s)
  for j=1:length(A)
    tau_1 = A(j,i);
    tau_2 = A(j,k);
    A(j,i) = c*tau_1 - s*tau_2;
    A(j,k) = s*tau_1 + c*tau_2;
  end
end