function multipliers = multipliers_pqz_schur(J)

% calculate multipliers
global lds

q = size(J,1)-1;

p1 = speye(q);
r = lds.col_coords;




J = J(1:q,1:q);

ntst   = lds.ntst;
ncol   = lds.ncol;
nphase = lds.nphase;
A = zeros(nphase, nphase, ntst);
B = zeros(nphase, nphase, ntst);
p_row_indices   = (1:nphase) + (ncol-1)*nphase;
j_row_indices   =  1:(ncol*nphase);
j_col_indices_p = (1:(ncol*nphase)) + nphase;
j_col_indices_A = (1:nphase);
j_col_indices_B = (1:nphase) + ncol*nphase;
for i=1:ntst
  sJ1 = J(r,lds.nphase+r);
  [sl1,~] = lu(sJ1);
  p1(r,r) = inv(sl1);
  r = r+lds.ncol_coord;
  sJ = J(j_row_indices,j_col_indices_p);
  [sl,~] = lu(sJ);
  p = inv(sl);
  p = p(p_row_indices,:);
  A(:,:,i) = p * J(j_row_indices,j_col_indices_A);
  B(:,:,i) = p * J(j_row_indices,j_col_indices_B);
  j_row_indices    = j_row_indices   + ncol * nphase;
  j_col_indices_p  = j_col_indices_p + ncol * nphase;
  j_col_indices_A  = j_col_indices_A + ncol * nphase;
  j_col_indices_B  = j_col_indices_B + ncol * nphase;
end
Ac = num2cell(A,[1 2]);
Bc = num2cell(B,[1 2]);
S = [blkdiag(Ac{:}) zeros(ntst * nphase, nphase)] + ...
    [zeros(ntst * nphase, nphase) blkdiag(Bc{:})];
S1 = p1(lds.multi_r1,:)*J(:,lds.multi_r2);


[A,B] = pqzschur(A,B);
multipliers = ones(lds.nphase,1);
for i=1:lds.ntst
  multipliers = multipliers .* diag(A(:,:,i)) ./ diag(B(:,:,i));
end
multipliers = sort(multipliers,'descend', 'ComparisonMethod', 'abs');  
print_diag(3,'multiplier with norm larger than 0.7: %.15f\n', ...
  multipliers(abs(multipliers)>0.7));

