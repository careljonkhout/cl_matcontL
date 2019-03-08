function multipliers = multipliers_explicit_p(J)

% calculate multipliers
global lds

q = size(J,1)-1;
J = J(1:q,1:q);
p = speye(q);
r = lds.col_coords;
number_of_nonzeros_p = lds.ntst*numel(lds.col_coords)^2;
p_row_indices = zeros(number_of_nonzeros_p,1);
p_col_indices = zeros(number_of_nonzeros_p,1);
p_values = zeros(number_of_nonzeros_p,1);
p_element_counter = 0;
for i=lds.tsts
  sJ = J(r,lds.nphase+r);
  [sl,~] = lu(sJ);
  m = inv(sl);
  for j=r
    range = p_element_counter + (1:numel(r));
    p_row_indices(range) = j * ones(numel(r),1);
    p_col_indices(range) = r;
    p_values(range) = m(j-r(1)+1,:);
    p_element_counter = p_element_counter + numel(r);
  end
  r = r+lds.ncol_coord;
end
p = sparse(p_row_indices,p_col_indices,p_values,q,q);
S = p(lds.multi_r1,:)*J(:,lds.multi_r2);

mon_straightforward = eye(lds.nphase);

for i=lds.ntst:-1:1
  indices1 = (1:lds.nphase) + (i-1)*lds.nphase;
  indices2 = indices1 + lds.nphase;
  mon_straightforward = mon_straightforward ...
    / S(indices1,indices2) * S(indices1,indices1); 
end



multipliers = eig(mon_straightforward);
multipliers = sort(multipliers,'descend', 'ComparisonMethod', 'abs'); 



