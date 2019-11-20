function [file] = get_latest_singularity_file(path, base_filename)
  if nargin < 2
    error(['get_latest_singularity_file needs two arguments: ' ...
           '''path'' and ''basename'''])
  end
  listing       = dir(fullfile(path, 'Data', [base_filename '*.mat']));
  if isempty(listing)
    error('get_latest_singularity_file:not_found', ...
          'No .mat files were found in in the Data folder starting with %s', ...
                                        base_filename);
  end
  listing_table = struct2table(listing);
  listing_table = sortrows(listing_table, 'name');
  listing       = table2struct(listing_table);
  filename      = listing(end).name;
  file          = fullfile(path, 'Data', filename);