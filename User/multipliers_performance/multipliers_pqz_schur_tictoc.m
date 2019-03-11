function multipliers = multipliers_pqz_schur_tictoc(J)

% calculate multipliers
global lds
tic
q = size(J,1)-1;
J = J(1:q,1:q);
p = speye(q);
r = lds.col_coords;
for i=lds.tsts
  sJ = J(r,lds.nphase+r);
  [sl,~] = lu(sJ);
  p(r,r) = inv(sl);
  r = r+lds.ncol_coord;
end
S = p(lds.multi_r1,:)*J(:,lds.multi_r2);
toc
A = zeros(lds.nphase,lds.nphase,lds.ntst);
B = zeros(lds.nphase,lds.nphase,lds.ntst);

for i=1:lds.ntst
  indices1 = (1:lds.nphase) + (i-1)*lds.nphase;
  indices2 = indices1 + lds.nphase;
  A(:,:,i) = S(indices1,indices1);
  B(:,:,i) = S(indices1,indices2);
end
[A,B] = pqzschur(A,B);
multipliers = ones(lds.nphase,1);
for i=1:lds.ntst
  multipliers = multipliers .* diag(A(:,:,i)) ./ diag(B(:,:,i));
end
multipliers = sort(multipliers,'descend', 'ComparisonMethod', 'abs');  
print_diag(3,'multiplier with norm larger than 0.7: %.15f\n', ...
  multipliers(abs(multipliers)>0.7));

