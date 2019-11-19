% Print something to the Matlab console. The second time (and any subsequent 
% times) print_temp is called it removes the text printed the previous time by
% printing n_bytes backspace characters, before printing the new text.
function print_temp(formatstring, varargin)
  persistent n_bytes
  if ~ isempty(n_bytes)
    fprintf(repmat('\b', 1, n_bytes))
  end
  if nargin > 0
    n_bytes = fprintf(formatstring, varargin{2:end});
  else
    n_bytes = [];
  end
end
  