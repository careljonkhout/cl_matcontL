function V = NP_SS_compute_subspace(period, parameters)
  global cds
  
  cds.mv_count = 0;
  
  monodromy_map = @(x) NP_SS_monodromy_map(x, period, parameters);
  
  using_octave_eigs = true;
  
  if using_octave_eigs
    eigensolver_limit = cds.n_phases / 2;
    nEigs = min(eigensolver_limit, cds.p + 1);
    options.tol = 1e-12;
    [eigenvectors, eigenvalue_matrix] = eigs(monodromy_map, cds.n_phases, ...
                                           nEigs, 'LM', options);
                                           
  
    eigenvalues = diag(eigenvalue_matrix);
    [eigenvalues, sorting_permutation] = sort(eigenvalues, 'descend');
    eigenvectors = eigenvectors(:, sorting_permutation);
      print_diag(2,'computing subspace mv_count: %d\n', cds.mv_count);
  
    eigenvalues = diag(eigenvalues);
    basis = zeros(cds.n_phases, cds.p + 1);

    i = 0;

    while i <= cds.p && i <= eigensolver_limit - 1
      i = i + 1;
      basis(:, i) = real(eigenvectors(:, i));
      if abs(imag(eigenvalues(i))) > 0
        i = i + 1;
        basis(:, i) = imag(eigenvectors(:, i-1));
      end
    end
    V = orth(basis(:,1:i));
  end
  
  using_ahbschur = false;
  
  
  
  
  if using_ahbschur
    n_eigenvalues = min(cds.n_phases, cds.preferred_basis_size + 2);
    options.k = 15;
    options.blsz = 1;

    options.tol = 1e-30;
    [Q, T, convergence] = ahbschur(monodromy_map, ...
                                                 cds.n_phases, ...
                                                 options);
    
                    
    basis = zeros(cds.n_phases, cds.preferred_basis_size + 1);

    j = 1;

    while j < 14
      basis(:, j) = Q(:, j);
      j = j + 1;
      if T(j, j-1) ~= 0
        basis(:, j) = Q(:, j + 1);
        j = j + 1;
      end
    end
    V = basis(:,1:j-1);                               
                                                 
  end
                                                 
