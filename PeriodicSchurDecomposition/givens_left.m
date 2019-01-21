%%
% Compute the matrix G * A where G is a givens rotation
% the following would give the same result (but is less efficient):
% givens_mat = eye(size(A));
% givens_mat(i,i) = c;
% givens_mat(k,k) = c;
% givens_mat(i,k) = s;
% givens_mat(k,i) = -s;
% A = givens_mat*A
%
function A = givens_left(A,i,k,c,s)
  for j=1:length(A)
    tau_1 = A(i,j);
    tau_2 = A(k,j);
    A(i,j) =   c*tau_1 + s*tau_2;
    A(k,j) = - s*tau_1 + c*tau_2;
  end
end