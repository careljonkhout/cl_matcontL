function multipliers = multipliers(J)

% calculate multipliers
global lds

ntst   = lds.ntst;
ncol   = lds.ncol;
nphase = lds.nphase;

J_blocks = zeros(nphase*ncol,nphase*(ncol+1),ntst);
A        = zeros(nphase,     nphase,         ntst);
B        = zeros(nphase,     nphase,         ntst);

% retrieve block of jacobian first and store them as dense matrices
% to prevent repeating expensive lookups of elements in sparse matririces
j_row_indices = 1:ncol*nphase;
j_col_indices = 1:((ncol+1)*nphase);
for i=1:ntst
  J_blocks(:,:,i) = J(j_row_indices,j_col_indices);
  j_row_indices = j_row_indices + ncol*nphase;
  j_col_indices = j_col_indices + ncol*nphase;
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


[A,B] = pqzschur(A,B);
multipliers = ones(lds.nphase,1);
for i=1:lds.ntst
  multipliers = multipliers .* diag(A(:,:,i)) ./ diag(B(:,:,i));
end
multipliers = sort(multipliers,'descend', 'ComparisonMethod', 'abs');

critical_multipliers = multipliers(abs(multipliers)>0.85);
text = 'multipliers with norm larger than 0.85:\n';
for i=1:length(critical_multipliers)
  m = critical_multipliers(i);
  if abs(imag(m)) > 1e-8
    % Call me crazy, but I'd like a space between a + or - and a number.
    % This is the only way I know how to do this.
    if  sign(imag(m)) > 1
      sign_text = '+';
    else
      sign_text = '-';
    end
    % This text will not be so long that it will negatively impact performance.
    % Therefore, the 'array grows on every loop iteration'-warning is ignored
    text  = [ text sprintf('%.15f %s %.15fi\n', ...
              real(m), sign_text, abs(imag(m))) ]; %#ok<AGROW>
  else
    text  = [ text sprintf('%.15f\n', real(m)) ]; %#ok<AGROW>
  end
end

print_diag(3,text);

