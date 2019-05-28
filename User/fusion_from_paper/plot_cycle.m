load('/home/carel/Documents/cl_matcontL/User/fusion_from_paper/Selected_Data/fusion_from_paper_oc_2019_05_25_17_29_4.595/point_00000001.mat')
% number of "phase variables" (dependent variables of the system of ODEs):
nphase = numel(point.multipliers); 
ntst = point.ntst; % number of mesh intervals
ncol = point.ncol; % number of collocation points
fine_mesh = get_fine_mesh(point.timemesh, ntst, ncol);
x = reshape(point.x(1:end-2), nphase, ntst * ncol +1)';
coords_to_plot = 3;
x = x(:, coords_to_plot);
index = 64*4;
rotated_time = fine_mesh(index:end) - fine_mesh(index);
rotated_time = [rotated_time (rotated_time(end) + fine_mesh(1:index))];
rotated_x    = x(index:end);
rotated_x    = [rotated_x; x(1:index)];

%rotated_time = fine_mesh;
%rotated_x    = x;


plot(rotated_time*point.T, rotated_x);
xlim([0 point.T]);
xlabel('time')
ylabel('Z_0')