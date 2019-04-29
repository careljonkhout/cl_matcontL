function file = get_latest_mat_file(dirname)
  listing       = dir(fullfile(dirname,'*.mat'));
  listing_table = struct2table(listing);
  listing_table = sortrows(listing_table, 'datenum');
  listing       = table2struct(listing_table);
  file          = listing(end).name;