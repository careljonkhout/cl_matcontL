load('/home/carel/Documents/cl_matcontL/User/fusion_from_paper/Data/fusion_from_paper_oc_03-May-2019_13_13_30/point_00000011.mat')

ntst = point.ntst;
ncol = point.ncol;

fine_mesh = zeros(1,ntst*ncol+1);
fine_mesh_index = 2;

for i=1:ntst
  mesh_interval_width = point.timemesh(i+1) - point.timemesh(i);
  fine_mesh_dt        = mesh_interval_width / ncol;
  for j=1:ncol
    fine_mesh(fine_mesh_index) = point.timemesh(i) + j * fine_mesh_dt;
    fine_mesh_index            = fine_mesh_index + 1;
  end
end

period = point.x(end-1);

coordinate_data = reshape(point.x(1:end-2),147,length(fine_mesh));

t = linspace(0, period, 1000);

y = interp1(period*fine_mesh,(coordinate_data - coordinate_data(:,1))', t, 'spline');

plot(t, y)