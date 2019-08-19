function result = emat2(template_filename)
  template = fileread(template_filename);
  
  block_starts = strfind(template, '<%');
  block_ends   = strfind(template, '%>');
  
  full_blocks  = cell(length(block_starts),1);
  eval_results = cell(length(block_starts),1);
  
  for i = 1:length(block_starts)
    full_blocks{i}  = template(block_starts(i)  :block_ends(i)+1);
    block           = template(block_starts(i)+3:block_ends(i)-2);
    eval_results{i} = to_string(evalin('caller', block));
  end
  
  result = replace_symbols(template, full_blocks, eval_results);
end

function str = to_string(obj)
  if ischar(obj)
    str = obj;
  elseif isnumeric(obj) && floor(obj) == obj
    str = num2str(obj);
    error('block in template returned unexpected result');
  end
end