function multipliers = multipliers_old(J)

% calculate multipliers
global lds
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
S = full(p(lds.multi_r1,:)*J(:,lds.multi_r2));


A = zeros(lds.nphase,lds.nphase,lds.ntst);
B = zeros(lds.nphase,lds.nphase,lds.ntst);

for i=1:lds.ntst
  indices1 = (1:lds.nphase) + (i-1)*lds.nphase;
  indices2 = indices1 + lds.nphase;
  A(:,:,i) = S(indices1,indices1);
  B(:,:,i) = S(indices1,indices2);
  
end
[A,B] = pqzschur(A,B);
eig_pqzschur = ones(lds.nphase,1);
for i=1:lds.ntst
  eig_pqzschur = eig_pqzschur .* diag(A(:,:,i)) ./ diag(B(:,:,i));
end
eig_pqzschur = sort(eig_pqzschur,'descend', 'ComparisonMethod', 'abs');  
print_diag(3,'eig_pqzschur: %.15f\n', eig_pqzschur(1:5));



for i=(1:(lds.ntst-1))*lds.nphase
  r = i+lds.phases;
  for j=r
    f = S(r,j)/S(j-lds.nphase,j);
    S(r,:) = S(r,:)-f*S(j-lds.nphase,:);
  end
end
r1 = (lds.ntst-1)*lds.nphase+lds.phases;
A0 = S(r1,lds.phases);
A1 = S(r1,r1+lds.nphase);


multipliers = eig(-A0,A1);
multipliers = sort(multipliers,'descend', 'ComparisonMethod', 'abs');

lds.monodromy = -A1\A0;

multipliers_mon = eig(lds.monodromy);
multipliers_mon = sort(multipliers_mon,'descend', 'ComparisonMethod', 'abs');
print_diag(3,'multiplier_mon: %.15f\n', multipliers_mon(1:5));

