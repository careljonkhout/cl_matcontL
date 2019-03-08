function multipliers = multipliers(J)

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
S = p(lds.multi_r1,:)*J(:,lds.multi_r2);


% Carel: In case I forget to remove it, code form here until line 31 can be
% deleted

mon_straightforward = eye(lds.nphase);

for i=lds.ntst:-1:1
  indices1 = (1:lds.nphase) + (i-1)*lds.nphase;
  indices2 = indices1 + lds.nphase;
  mon_straightforward = mon_straightforward ...
    / S(indices1,indices2) * S(indices1,indices1); 
end



multipliers = eig(mon_straightforward);
multipliers = sort(multipliers,'descend', 'ComparisonMethod', 'abs'); 



