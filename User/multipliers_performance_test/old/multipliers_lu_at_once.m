function multipliers = multipliers(J)

% calculate multipliers
global lds contopts
ntst = lds.ntst;
ncol = lds.ncol;
nphase = lds.nphase;

disp('gives wrong result');
q = size(J,1)-1;
J = J(1:q,1:q);
arg_for_lu = J;
arg_for_lu(lds.mask==0) = 0;
arg_for_lu = circshift(arg_for_lu,lds.nphase,2);
[l,~]=lu(arg_for_lu);
p = inv(l);
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

