function out_str=replace_symbols(in_str, old_symbols, new_symbols)
  parselist = parse_expression(in_str, old_symbols);
  out_str = '';
  for i=1:length(parselist)
    if parselist{i}.is_symbol
      out_str = ...
        [out_str new_symbols{parselist{i}.symbol_index}]; %#ok<AGROW>
    else
      out_str = [out_str parselist{i}.data]; %#ok<AGROW>
    end
  end       
end

% str must be a char array
% symbols must be a cell array of horizontal char arrays
function parse_list = parse_expression(str, symbols)
  parse_list = {};
  parse_list{1}.is_symbol = false;
  parse_list{1}.data = str;
  for si = 1:length(symbols) % si means symbols_index
    pli = 1; % i.e. parse_list_index
    while pli <= length(parse_list) 
      if parse_list{pli}.is_symbol
        pli = pli + 1;
      else
        parsed_sub_string = split_at_symbol( ...
            parse_list{pli}.data, symbols{si}, si);
        parse_list = {parse_list{1:pli-1} ...
            parsed_sub_string{1:end} parse_list{pli+1:end}};
        pli = pli + length(parsed_sub_string);
      end
    end
  end
end

function parse_list = split_at_symbol(str, symbol, symbol_index)
  parse_list = {};
  symbol_indices = strfind(str, symbol);     
  % select only matches that are isolated
  % i.e. not part of a larger symbol:
  isolated = arrayfun( @(si) ...
      (si == 1 ...
    || ~ isletter(str(si-1))) ...
    && (si + length(symbol) > length(str) ...
    ||  ~ isstrprop(str(si+length(symbol)), 'alphanum')), ...
        symbol_indices);

  symbol_indices = symbol_indices(isolated);
  % if symbol is not found in str
  % return a trivial parse_list
  if isempty(symbol_indices)
    parse_list{1}.is_symbol = false;
    parse_list{1}.data = str;
    return;
  % if str does not start with symbol
  % put the part of the string up to the first symbol
  % in the first node of the parselist
  elseif symbol_indices(1) > 1
    parse_list{1}.is_symbol = false;
    parse_list{1}.data = str(1:symbol_indices(1)-1);
  end
  for i=1:length(symbol_indices)
    parse_list{end+1}.is_symbol = true; %#ok<AGROW>
    parse_list{end}.data = ...
        str(symbol_indices(i):symbol_indices(i)+length(symbol)-1);
    parse_list{end}.symbol_index = symbol_index;
    % we process the part of str after the current location of symbol
    if i < length(symbol_indices)
        parse_list{end+1}.is_symbol = false; %#ok<AGROW>
        parse_list{end}.data = str( ...
          symbol_indices(i)+length(symbol):symbol_indices(i+1)-1);
    % after the last occurence of symbol is processed
    % we put the remainder of str in the last node of parse_list
    else
      parse_list{end+1}.is_symbol = false; %#ok<AGROW>
      parse_list{end}.data = str(symbol_indices(i)+length(symbol):end);            
    end
  end
end

% for debugging
% one might also use celldisp for prettier output
function print_parselist(parse_list)
  for e_cell=parse_list
    e = e_cell{1};
    if e.is_symbol
      fprintf("%s %d",e.data,e.symbol_index);
    else
      fprintf("'%s'",e.data);
    end
  end
end
