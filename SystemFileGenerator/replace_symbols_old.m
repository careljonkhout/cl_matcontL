function out_str=replace_symbols(in_str, old_symbols, new_symbols)
  parse_list = parse_expression(in_str, old_symbols);
  out_str_parts = cell(length(parse_list),1);
  for i=1:length(parse_list)
    if parse_list{i}.is_symbol
      out_str_parts{i} = new_symbols{parse_list{i}.symbol_index};
    else
      out_str_parts{i} = parse_list{i}.data;
    end
  end
  out_str = strjoin(out_str_parts, '');
end

% str must be a char array
% symbols must be a cell array of horizontal char arrays
function parse_list = parse_expression(str, symbols)
  parse_list = {};
  parse_list{1}.is_symbol = false;
  parse_list{1}.data = str;
  for si = 1:length(symbols) % si means symbols_index
    next_parse_lists = cell(length(parse_list),1);
    for pli = 1:length(parse_list) % i.e. parse_list_index
      if parse_list{pli}.is_symbol
        next_parse_lists(pli) = {parse_list(pli)};
      else
        parsed_sub_string = split_at_symbol( ...
            parse_list{pli}.data, symbols{si}, si);
        next_parse_lists{pli} = parsed_sub_string;
      end
    end
    parse_list = flatten(next_parse_lists);
  end
end

function parse_list = split_at_symbol(str, symbol, symbol_index)
  parse_list = {};
  match_indices = strfind(str, symbol);     
  % select only matches that are isolated
  % i.e. not part of a larger symbol:
  isolated = arrayfun( @(mi) is_isolated(str, symbol, mi), match_indices);
  match_indices = match_indices(isolated);
  % if symbol is not found in str
  % return a trivial parse_list
  if isempty(match_indices)
    parse_list{1}.is_symbol = false;
    parse_list{1}.data = str;
    return;
  % if str does not start with symbol
  % put the part of the string up to the first symbol
  % in the first node of the parse_list
  elseif match_indices(1) > 1
    parse_list{1}.is_symbol = false;
    parse_list{1}.data = str(1:match_indices(1)-1);
  end
  for i=1:length(match_indices)
    parse_list{end+1}.is_symbol = true; %#ok<AGROW>
    parse_list{end}.data = ...
        str(match_indices(i):match_indices(i)+length(symbol)-1);
    parse_list{end}.symbol_index = symbol_index;
    % we process the part of str after the current location of symbol
    if i < length(match_indices)
        parse_list{end+1}.is_symbol = false; %#ok<AGROW>
        parse_list{end}.data = str( ...
          match_indices(i)+length(symbol):match_indices(i+1)-1);
    % after the last occurence of symbol is processed
    % we put the remainder of str in the last node of parse_list
    elseif match_indices(i) + length(symbol) <= length(str)
      parse_list{end+1}.is_symbol = false; %#ok<AGROW>
      parse_list{end}.data = str(match_indices(i)+length(symbol):end);            
    end
  end
end

function is_isolated = is_isolated(str, symbol, match_index)
  if (match_index == 1)
    is_isolated_left = true;
  else
    prev_char = str(match_index - 1);
    is_isolated_left = ~ isletter(prev_char) && prev_char ~= '_';
  end
 
  if match_index + length(symbol) == length(str) + 1
    is_isolated_right = true;
  else
    next_char = str(match_index + length(symbol));
    is_isolated_right = ~ isstrprop(next_char, 'alphanum') && next_char ~= '_';
  end
  is_isolated = is_isolated_left && is_isolated_right;
end

% for debugging
% one might also use celldisp for prettier output
function print_parse_list(parse_list) %#ok<DEFNU>
  for e_cell=parse_list
    e = e_cell{1};
    if e.is_symbol
      fprintf("%s %d",e.data,e.symbol_index);
    else
      fprintf("'%s'",e.data);
    end
  end
end

function flattened = flatten(array)
  new_array_length = 0;
  for i=1:length(array)
    new_array_length = new_array_length + length(array{i});
  end
  flattened = cell(new_array_length,1);
  flattened_index = 1;
  for i=1:length(array)
    flattened(flattened_index:flattened_index + length(array{i})-1) = array{i};
    flattened_index = flattened_index + length(array{i});
  end
end