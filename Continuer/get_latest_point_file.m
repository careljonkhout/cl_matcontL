function [file, point_index] = get_latest_point_file(dirname)
  listing       = dir(fullfile(dirname,'*.mat'));
  if isempty(listing)
    error('No .mat files were found in %s.', dirname)
  end
  listing_table = struct2table(listing);
  listing_table = sortrows(listing_table, 'datenum');
  listing       = table2struct(listing_table);
  filename      = listing(end).name;
  point_index   = sscanf(filename, 'point_%d');
  if isempty(point_index)
    error('No file whose name starts with point_ was found in %s', dirname);
  end
  file          = fullfile(dirname, filename);