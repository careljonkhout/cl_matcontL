function for_each_two_consecutive_points_in(dirname, func)
  if ~ isa(func, 'function_handle')
    error(['second arugument of ''for_each_point_in'' should be ' ...
            'a function handle']);
  end
  files       = dir(fullfile(dirname,'point_*.mat'));
  if isempty(files)
    error('No point_*.mat files were found in %s.', dirname)
  end
  listing_table = struct2table(files);
  listing_table = sortrows(listing_table, 'name');
  files         = table2struct(listing_table);
  point_2       = load_point(files, 1);
  for i = 2 : length(files)
    point_1 = point_2;
    point_2 = load_point(files, i);
    func(point_1, point_2);
  end
end

function point = load_point(files, i)
  load(fullfile(files(i).folder, files(i).name), 'point');
end
  