N               = 75;
fprintf('N=%d\n',N);
odefile         = str2func(sprintf('fusion_precomputed_with_sage_N_%d',N));
syms              a b q_inf
xxxxx           = sym('xxxxx', [1 3*(N-1)]);
handles             = feval(odefile);
jacobian_function   = handles{3};
jacobian = feval(jacobian_function,0,xxxxx,a,b,q_inf);
jacobian_c_code     = ccode(jacobian);
nphases = 3*(N-1);
for i=nphases:-1:1
  find    = sprintf('xxxxx%d',i);
  replace = sprintf('xxxxx[%d]',i-1);
  jacobian_c_code = strrep(jacobian_c_code, find, replace);
end
lines   = strsplit(jacobian_c_code, '\n');
pattern = 'jacobian\[(\d+)\]\[(\d+)\] = (.*);';
% array to store structs with data for jacobian elements
e = cell(length(lines),1); 
for i=1:length(lines)
  tokens     = regexp(lines{i}, pattern, 'tokens', 'all');
  e{i}.row_index = str2double(tokens{1}{1});
  e{i}.col_index = str2double(tokens{1}{2});
  e{i}.expression = tokens{1}{3};
  e{i}.linear_index = e{i}.col_index * nphases + e{i}.row_index;
end


linear_indices = cellfun(@(x) x.linear_index, e);
% find the permutation that sorts the jacobian elements by column index
[~, permutation] = sort(linear_indices);
% apply the permutation  to sort the jacobian elements by column index
e = e(permutation);

% caution: the next while loop assumes that there are no empty columns

jc = zeros(nphases,1);
c_style_e_index = 0;
for c_style_col_index = 0:nphases-1  
  jc(c_style_col_index+1) = c_style_e_index;
  while c_style_e_index+1 <= length(lines) && ...
      e{c_style_e_index+1}.col_index == c_style_col_index
    c_style_e_index = c_style_e_index + 1;
  end
end
jc(nphases+1) = length(lines);
jc_str = '{';

for i=1:length(jc)
  jc_str = [jc_str sprintf('%d, ', jc(i)) ]; %#ok<AGROW>
end
jc_str(end+1) = '}';

filename        = sprintf('fusion_sparse_jacobian_c_code_N_%d',N);
fid             = fopen(filename,'w');
%fprintf(fid, 'jc = {');
%fprintf(fid, '%d, ', jc);
%fprintf(fid, '};\n\n');

%fprintf(fid, 'row_indices = {');
%fprintf(fid, '%d, ', cellfun(@(x) x.row_index, e));
%fprintf(fid, '};\n\n');

expressions = cellfun(@(x) x.expression, e, 'UniformOutput', false);

jacobian_values = '';
row_indices = '{';

for i=1:length(lines)
  jacobian_values = [jacobian_values ...
      sprintf('jacobian_values[%d] = %s;\n', i-1, e{i}.expression)]; %#ok<AGROW>
  row_indices = [row_indices sprintf('%d, ',      e{i}.row_index )]; %#ok<AGROW>
end
row_indices(end+1) = '}';

nzmax= length(lines);


c_code = emat('fusion_sparse_jacobian_template.c.emat');

fprintf(fid,c_code);


fclose(fid);

[status, cmd_out] = ...
  system(sprintf('sh ./compile_mex_fusion_sparse_jac_%d.sh',N));

if status ~= 0
  disp(cmd_out)
end

test_fusion_mex_quick