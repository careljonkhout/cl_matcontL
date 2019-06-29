function points = load_matfile_points
  files = dir;
  for i = 1:length(files)
    file = files(i);
    if length(file.name) >= 5 && strcmp(file.name(1:5),'point')
      point_index = sscanf(file.name, 'point_%d.mat');
      load(file.name, 'point');
      points{point_index} = point; %#ok<AGROW>
                       % preallocating does not appear to be worth the effort
    end
  end
  if ~ exist('points', 'var')
    error('No points saved as matfiles were found in the current folder.');
  end
end
      
      
  