function multipliers = multipliers_pqz_schur(J)

% calculate multipliers
global lds

q = size(J,1)-1;
J = J(1:q,1:q);
A = zeros(lds.nphase,lds.nphase,lds.ntst);
B = zeros(lds.nphase,lds.nphase,lds.ntst);
ncol   = lds.ncol;
nphase = lds.nphase;
p_row_indices   = (1:nphase) + (ncol-1)*nphase;
j_row_indices   =  1:(ncol*nphase);
j_col_indices   = (1:nphase);
for i=1:lds.ntst
  sJ = J(j_row_indices,j_row_indices + nphase);
  [sl,~] = lu(sJ);
  p = inv(sl);
  p = p(p_row_indices,:);
  A(:,:,i) = p * J(j_row_indices,j_col_indices);
  B(:,:,i) = p * J(j_row_indices,j_col_indices + ncol*nphase);
  j_col_indices  = j_col_indices + ncol * nphase;
  j_row_indices  = j_row_indices + ncol * nphase;
end


[A,B] = pqzschur(A,B);
multipliers = ones(lds.nphase,1);
for i=1:lds.ntst
  multipliers = multipliers .* diag(A(:,:,i)) ./ diag(B(:,:,i));
end
multipliers = sort(multipliers,'descend', 'ComparisonMethod', 'abs');  
print_diag(3,'multiplier with norm larger than 0.7: %.15f\n', ...
  multipliers(abs(multipliers)>0.7));

