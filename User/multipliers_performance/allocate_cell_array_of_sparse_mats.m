function result = allocate_cell_array_of_sparse_mats(m,n,p)
  result = cell(p,1);
  for i=1:p
    result{i} = spalloc(m,n,m);
  end