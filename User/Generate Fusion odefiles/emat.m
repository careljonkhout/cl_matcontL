function result = emat(template_file)
  template = fileread(template_file);
  
  block_starts = strfind(template, '<%');
  block_ends = strfind(template, '%>');
  
  full_blocks  = cell(length(block_starts),1);
  eval_results = cell(length(block_starts),1);
  
  for i = 1:length(block_starts)
    full_blocks{i} = template(block_starts(i):block_ends(i)+1);
    block = template(block_starts(i)+4:block_ends(i)-2);
    eval_results{i} = evalin('caller', block);
  end
  
  result = replace_symbols(template, full_blocks, eval_results);
end
