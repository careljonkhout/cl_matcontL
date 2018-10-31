function add_fusion_systems_to_path()
  filename = mfilename;
  fullpath = mfilename('fullpath');
  path_wo_filename = fullpath(1:end-length(filename));
  fusion_systems_path = fullfile(path_wo_filename, ...
    '..', '..', 'Systems', 'fusion');
  addpath(fusion_systems_path);