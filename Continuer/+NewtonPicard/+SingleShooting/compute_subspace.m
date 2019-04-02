function V = compute_subspace(period, parameters)
  global cds
  
  p = min([cds.preferred_basis_size cds.nphases]);
  cds.p_extra = 2;
  cds.p = p;
  cds.mv_count = 0;
  
   monodromy_map = @(x) NewtonPicard.SingleShooting.monodromy_map( ...
                        x, period, parameters);
  
  [eigenvectors, eigenvalues, no_convergence] = eigs(monodromy_map, ...
                                                  cds.nphases, p + cds.p_extra);
  print_diag(0,'computing subspace mv_count: %d\n', cds.mv_count);
  if no_convergence
    V = [];
    fprintf(['Newton_Picard_Correction.m:', ...
      ' eigenvalues of monodromy matrix did not converge.\n'])
    return
  end

  eigenvalues = diag(eigenvalues);
  basis = zeros(cds.nphases, p + cds.p_extra);
  cds.eigenvalues = eigenvalues;
  i = 0;

  while i <= p + cds.p_extra - 1
    i = i + 1;
    basis(:, i) = real(eigenvectors(:,i));
    if abs(imag(eigenvalues(i))) > 0
      i = i + 1;
      basis(:, i) = imag(eigenvectors(:,i-1));
    end
  end
  V = orth(basis(:,1:i));