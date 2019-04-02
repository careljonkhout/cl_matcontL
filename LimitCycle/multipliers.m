function multipliers = multipliers(J)

% calculate multipliers
global lds contopts

ntst   = lds.ntst;
ncol   = lds.ncol;
nphase = lds.nphase;

J_blocks = zeros(nphase*ncol,nphase*(ncol+1),ntst);
A        = zeros(nphase,     nphase,         ntst);
B        = zeros(nphase,     nphase,         ntst);

% retrieve block of jacobian first and store them as dense matrices
% to prevent repeating expensive lookups of blocks in sparse matririces
j_row_indices = 1 :  ncol    * nphase;
j_col_indices = 1 : (ncol+1) * nphase;
for i=1:ntst
  J_blocks(:,:,i) = J(j_row_indices, j_col_indices);
  j_row_indices = j_row_indices + ncol * nphase;
  j_col_indices = j_col_indices + ncol * nphase;
end


p_row_indices   = (1:nphase) + (ncol-1)*nphase;
j_col_indices_p = (1:(ncol*nphase)) + nphase;
j_col_indices_A = (1:nphase);
j_col_indices_B = (1:nphase) + ncol*nphase;
for i=1:ntst
  sJ = J_blocks(:,j_col_indices_p,i);
  [sl,~] = lu(sJ);
  p = inv(sl);
  p = p(p_row_indices,:);
  
  A(:,:,i) = p * J_blocks(:,j_col_indices_A,i);
  B(:,:,i) = p * J_blocks(:,j_col_indices_B,i);
end

try
  [A,B] = pqzschur(A,B);
  multipliers = ones(lds.nphase,1);
  for i=1:lds.ntst
    multipliers = multipliers .* diag(A(:,:,i)) ./ diag(B(:,:,i));
  end
catch
  print_diag(3, 'pqzschur did not converge, using Gaussian elimination\n');
 
  S = zeros(ntst * nphase, (ntst+1) * nphase + 1);
  for i=1:ntst
    row_indices   = (1:nphase) + (i-1) * nphase;
    col_indices_A = (1:nphase) + (i-1) * nphase;
    col_indices_B = (1:nphase) +  i    * nphase;
    S(row_indices, col_indices_A) = A(:,:,i);
    S(row_indices, col_indices_B) = B(:,:,i);
  end
  
  for i = (1:ntst-1) * nphase
    r = i + (1:nphase);
    for j = r
      f      = S(r,j) /     S(j - nphase, j);
      S(r,:) = S(r,:) - f * S(j - nphase, :);
    end
  end
  row_indices    = (ntst-1) * nphase + (1:nphase);
  col_indices_A0 =                     (1:nphase);
  col_indices_A1 =  ntst    * nphase + (1:nphase);
  A0 = S(row_indices, col_indices_A0);
  A1 = S(row_indices, col_indices_A1);
  multipliers = eig(-A0,A1);
end

multipliers = sort(multipliers, 'descend', 'ComparisonMethod', 'abs');

critical_multipliers = ...
  multipliers(abs(multipliers) > contopts.multiplier_print_threshold);
text = sprintf('multipliers with norm larger than %.3f:\n', ...
                contopts.multiplier_print_threshold);
text = [text multipliers2str(critical_multipliers)];

print_diag(3,text);

