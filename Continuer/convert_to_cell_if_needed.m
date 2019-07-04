function out = convert_to_cell_if_needed(in)
  if ~ iscell(in)
    out = num2cell(in);
  else
    out = in;
  end
end