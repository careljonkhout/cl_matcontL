% Print the status of something that changes continuously to the Matlab console.
% The second time (and any subsequent times) print_temp is called it removes the
% text printed the previous time by printing n_bytes backspace characters,
% before printing the new text.

function print_temp(format_string, varargin)
  persistent n_bytes
  if ~ isempty(time_print_temp())
    dv = datevec(now - time_print_temp());
    % do not print anything if something was printed less than 0.1 seconds ago
    if dv(6) < 0.1
      return
    end
  end
  time_print_temp(now);
  % if something was printed before, remove it by printing backspaces ('\b')
  if ~ isempty(n_bytes)
    fprintf(repmat('\b', 1, n_bytes))
  end
  if nargin > 0
    n_bytes = fprintf(format_string, varargin{:});
  else
    n_bytes = [];
  end
end
  