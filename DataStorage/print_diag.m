function print_diag(priority, message, varargin)
%This function prints diagnostic messages  
global contopts cds

if ~ isnumeric(priority) || ~ any(ismember(priority, 0:6))
  error(['The first input ''priority'' to print_diag must be ' ...
         'an integer which is at least 0 and at most 6'])
end

if priority <= contopts.console_output_level
  fprintf(message,varargin{:});
end

if priority <= contopts.contL_DiagnosticsLevel
  if isfield(cds, 'logFID') && cds.logFID
    fprintf(cds.logFID, message,varargin{:});  
  end
end