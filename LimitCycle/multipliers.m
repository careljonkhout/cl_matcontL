function multipliers = multipliers(J)

% calculate multipliers
global lds

ntst   = lds.ntst;
ncol   = lds.ncol;
nphase = lds.nphase;

J_blocks = zeros(nphase*ncol,nphase*(ncol+1),ntst);
A        = zeros(nphase,     nphase,         ntst);
B        = zeros(nphase,     nphase,         ntst);

% Retrieve blocks of Jacobian first and store them as dense matrices to prevent
% repeating expensive lookups of blocks in sparse matririces.
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

% For some information on the methods used here see the article 
% ``improved numerical Floquet multipliers'' by Kurt Lust.


try
  [A,B] = pqzschur(A,B);
  multipliers = ones(lds.nphase,1);
  for i=1:lds.ntst
    multipliers = multipliers .* diag(A(:,:,i)) ./ diag(B(:,:,i));
  end
catch
  print_diag(3, 'pqzschur did not converge, using matrix multiplication\n');
  
  M = eye(nphase);
  for i=ntst:-1:1
    M = M  / B(:,:,i) * A(:,:,i);
  end
  multipliers = eig(M);
  
  if abs(multipliers(1)) > 10000
    % computing the multipliers by computing the monodromy matrix
    % straightforward matrix multication is not the most stable method. In the
    % article ``improved numerical Floquet multipliers'' by Kurt Lust, the
    % limitations of that approach are described. The method below, due to T.F.
    % Fairgrieve, described in the article 'OK Floquet Multipliers' is a bit
    % more stable, but this current implementation is slow, and uses a lot of
    % memory
    print_diag(3, 'large multipliers detected, using Gaussian elimination\n');

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
end

multipliers = sort(multipliers, 'descend', 'ComparisonMethod', 'abs');

print_diag(0,multipliers2str(multipliers));

