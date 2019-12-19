function [file, point_index] = get_latest_point_file(dirname)
  files       = dir(fullfile(dirname,'point_*.mat'));
  if isempty(files)
    error('No point_*.mat files were found in %s.', dirname)
  end
  listing_table = struct2table(files);
  listing_table = sortrows(listing_table, 'name');
  files         = table2struct(listing_table);
  filename      = files(end).name;
  point_index   = sscanf(filename, 'point_%d');
  file          = fullfile(dirname, filename);
end