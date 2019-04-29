


listing = dir(fullfile('Data','langford_ss','*.mat'));
listing_table = struct2table(listing);
listing_table = sortrows(listing_table, 'datenum');
listing = table2struct(listing_table);
struct2table(listing)
listing(end).name
