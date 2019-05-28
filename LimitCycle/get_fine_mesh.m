function fine_mesh = get_fine_mesh(coarse_mesh, ntst, ncol)
  fine_mesh = zeros(1, ntst * ncol + 1);
  fine_mesh_index = 2;

  for i=1:ntst
    mesh_interval_width = coarse_mesh(i+1) - coarse_mesh(i);
    fine_mesh_dt        = mesh_interval_width / ncol;
    for j=1:ncol
      fine_mesh(fine_mesh_index) = coarse_mesh(i) + j * fine_mesh_dt;
      fine_mesh_index            = fine_mesh_index + 1;
    end
  end
end