function print_diag(priority,message,varargin)
%This function prints diagnostic messages  
global contopts cds
DiagL = contopts.contL_DiagnosticsLevel;

if priority == 0 % always print priority 0 to the screen
%if priority >= 0 % always print priority 0 to the screen
    fprintf(message,varargin{:});
end

if priority <= DiagL
    if cds.logFID
        fprintf(cds.logFID,message,varargin{:});
    elseif priority ~= 0
        fprintf(message,varargin{:});
    end
end