function multipliers = multipliers(J)

% calculate multipliers
global lds contopts
ntst = lds.ntst;
ncol = lds.ncol;
nphase = lds.nphase;


q = size(J,1)-1;
J = J(1:q,1:q);
p = spalloc(q,q,lds.ntst*lds.nphase^2*lds.ncol);
r = lds.col_coords;
indices_inv_sl = (ncol-1)*nphase+1:ncol*nphase;
for i=lds.tsts
  sJ = J(r,lds.nphase+r);
  [sl,~] = lu(sJ);
  inv_sl = inv(sl);
  indices1 = (1:nphase) + (i-1)*nphase;
  p(indices1,r) = inv_sl(indices_inv_sl,:);
  r = r+lds.ncol_coord;
end
S = full(p*J(:,lds.multi_r2));
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

