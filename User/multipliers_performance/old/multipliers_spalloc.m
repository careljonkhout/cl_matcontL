function multipliers = multipliers(J)

% calculate multipliers
global lds contopts
q = size(J,1)-1;
J = J(1:q,1:q);
p = spalloc(q,q,lds.ntst*lds.nphase^2*lds.ncol^2);
r = lds.col_coords;
for i=lds.tsts
  sJ = J(r,lds.nphase+r);
  [sl,~] = lu(sJ);
  p(r,r) = inv(sl);
  r = r+lds.ncol_coord;
end
S = full(p(lds.multi_r1,:)*J(:,lds.multi_r2));
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

if contopts.nCriticalMultipliers > 0
  n = contopts.nCriticalMultipliers;
  multipliers = multipliers(1:n);
  % if we have a pair of complex multipliers, we want to include both. Hence if
  % the smallest multiplier in modulus of the set of critical multipliers (
  % which we define here ) is complex, and it's conjugate is not in the set of
  % critical multipliers, we disregard it by setting it to zero.
  if ~ isreal(multipliers(n)) && ...
      ~ is_a_conjugate_pair(multipliers(n), multipliers(n-1))
      multipliers(n) = 0;
  end
end
      

function is_a_conjugate_pair = is_a_conjugate_pair(a,b)
  is_a_conjugate_pair = abs(conj(a)-b) < 1e-14;
