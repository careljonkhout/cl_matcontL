
cd(get_path())
mex replace_symbols.c


argument_lists = { ...
    {}, ...        % replace_symbols requires 3 arguments, but 0 are given.
    {1}, ...       % replace_symbols requires 3 arguments, but 1 are given.
    {1 2}, ...     % replace_symbols requires 3 arguments, but 2 are given.
    {1 2 3 4}, ... % replace_symbols requires 3 arguments, but 2 are given.
    num2cell(ones(1000,1)), ...
    ...              % replace_symbols requires 3 arguments, but 1000 are given.
    {1 2 3}, ...        % The first argument is not a char array.
    {{} 2 3}, ...       % The first argument is not a char array.
    {int64(1) 2 3}, ... % The first argument is not a char array.
    {['a';'a'], 2, 3}, ... % 1st arg not one row
    {'a', 1, 1}, ...       % 2nd arg not cell array
    {'a', {'a'}, 1}, ...   % 3rd arg not cell array
    {'a', {'aa', 'aa'; 'bb', 'bb'}, {'a'}}, ... % 2nd arg not 1d
    {'a', {'a'}, {'aa', 'aa'; 'bb', 'bb'}}, ... % 3rd arg not 1d
    {'a', {1}, {'b'}}, ...   % non char in 2nd                 
    {'a', {'b'}, {1}}, ...   % non char in 3rd
  };
for i = 1:length(argument_lists)
  try 
    argument_list = argument_lists{i};
    replace_symbols(argument_list{:})
  catch my_error
    disp(my_error.message)
  end
end

disp(replace_symbols('a df',{'a'},{'b'}));