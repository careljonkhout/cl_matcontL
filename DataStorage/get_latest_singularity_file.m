function [file] = get_latest_singularity_file(dirname, base_file_name)
  listing       = dir(fullfile(dirname, 'Data', [base_file_name '*.mat']));
  if isempty(listing)
    error('No .mat files were found in in the Data folder starting with %s', ...
                                        base_file_name);
  end
  listing_table = struct2table(listing);
  listing_table = sortrows(listing_table, 'name');
  listing       = table2struct(listing_table);
  filename      = listing(end).name;
  file          = fullfile(dirname, 'Data', filename);