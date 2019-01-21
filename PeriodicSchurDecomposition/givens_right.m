%%
% Compute the matrix  A * G', where G is a givens rotation
% the following would give the same result (but is less efficient):
% givens_mat = eye(size(A));
% givens_mat(i,i) = c;
% givens_mat(k,k) = c;
% givens_mat(i,k) = s;
% givens_mat(k,i) = -s;
% A = A*givens_mat';
%
function A = givens_right(A,i,k,c,s)
  for j=1:length(A)
    tau_1 = A(j,i);
    tau_2 = A(j,k);
    A(j,i) =   c*tau_1 + s*tau_2;
    A(j,k) = - s*tau_1 + c*tau_2;
  end
end