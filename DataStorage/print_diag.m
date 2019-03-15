function print_diag(priority,message,varargin)
%This function prints diagnostic messages  
global contopts cds

if priority <= contopts.console_output_level
  fprintf(message,varargin{:});
end

if priority <= contopts.contL_DiagnosticsLevel && isfield(cds,'logFID') ...
  && cds.logFID
  fprintf(cds.logFID,message,varargin{:});  
end