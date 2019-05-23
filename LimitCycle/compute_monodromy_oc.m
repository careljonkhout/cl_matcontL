function monodromy = compute_monodromy_oc(x)
% calculate multipliers
global lds

ntst   = lds.ntst;
ncol   = lds.ncol;
nphase = lds.nphase;

A        = zeros(nphase,     nphase,         ntst);
B        = zeros(nphase,     nphase,         ntst);

p_row_indices   = (1:nphase) + (ncol-1)*nphase;
j_col_indices_p = (1:(ncol*nphase)) + nphase;
j_col_indices_A = (1:nphase);
j_col_indices_B = (1:nphase) + ncol*nphase;
for i=1:ntst
  block = collocation_block(x,i);
  sJ = block(:,j_col_indices_p);
  [sl,~] = lu(sJ);
  p = inv(sl);
  p = p(p_row_indices,:);
   
  A(:,:,i) = p * block(:,j_col_indices_A);
  B(:,:,i) = p * block(:,j_col_indices_B);
end

% For some information on the methods used here see the article 
% ``improved numerical Floquet multipliers'' by Kurt Lust.

monodromy = eye(nphase);
for i=ntst:-1:1 % from i = ntst downto 1
  monodromy = monodromy  / B(:,:,i) * A(:,:,i);
end
  