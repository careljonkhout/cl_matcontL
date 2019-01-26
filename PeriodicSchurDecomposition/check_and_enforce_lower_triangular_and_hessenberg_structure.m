function M = check_and_enforce_lower_triangular_and_hessenberg_structure(M)
  m = size(M,3);
  % check and enforce Hessenberg structure of M_m
  assert(all(all(abs(tril(M(:,:,m),-2))<eps(10))));
  M(:,:,m) = triu(M(:,:,m),-1);
  % check and enforce upper triangular structure of M_i for i=1:m-1
  for i=1:m-1
    assert(all(all(abs(tril(M(:,:,i),-1))<eps(10))));
    M(:,:,i) = triu(M(:,:,i));
  end
end

