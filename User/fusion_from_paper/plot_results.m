cd ~/Documents/cl_matcontL/User/fusion_from_paper
cd Selected_Data
cd fusion_from_paper_oc_2019_05_25_17_29_4.595
points = load_matfile_points;
points = points(1:300);
cd ..
load('fusion_from_paper_oc_2019_05_25_17_29_4.595.mat','s');

parameter_values = cellfun(@(c) c.x(end)  , points);
periods          = cellfun(@(c) c.x(end-1), points);

figure
hold on;
plot(parameter_values, periods);

singularities = s;

vertical_alignment = 'top';

for singularity = singularities
  if singularity.index <= 300
    plot_sing(singularity, 'VerticalAlignment', vertical_alignment)
    % switch vertical alignment after each singularity to reduce overlap of labels
    if strcmp(vertical_alignment, 'top')
      vertical_alignment = 'bottom';
    else
      vertical_alignment = 'top';
    end
  end
end


function plot_sing(s, varargin)
  L = s.data.x(end);
  T = s.data.x(end-1);
  plot(L,T,'r*')
  text(L,T,s.label,varargin{:})
end