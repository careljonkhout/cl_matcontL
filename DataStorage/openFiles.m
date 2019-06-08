function [datafile, logfile] = openFiles()
global cds contopts

% dbstack returns the function call stack. That is, denote the function that
% called this function by f_1, the function that called f_1 by f_2, the one that
% called f_2 by f_3, and so on. Then dbstack returns an array of structs with
% information about f_i in the i-th place in the array.
stack         = dbstack();
scriptname    = stack(end).file;
scriptname    = scriptname(1:end-2);
prob          = func2str(cds.probfile);
curve         = func2str(cds.curve);
% clock returns a numerical array in the following format:
% [<year> <month_number> <day> <hour> <minute> <seconds>}
% the seconds have five decimals after the comma
current_time  = clock;

if isempty(contopts.Filename) 
  %% create filename
  year    = current_time(1);
  month   = current_time(2);
  day     = current_time(3);
  hours   = current_time(4);
  minutes = current_time(5);
  seconds = current_time(6);
  strdate = sprintf( '%d_%02d_%02d',   year,  month,   day);
  % Milliseconds are nice to have in the filename, in case a run follows another
  % within one second.
  strtime = sprintf('%0d_%02d_%02.3f', hours, minutes, seconds);
  
  filename = [scriptname '_' strdate '_' strtime];
else
  filename = contopts.Filename;
end
cds.runID = filename;

%% Find path
if isempty(contopts.TestPath)   
  % note that it is attempted to set testpath when contset() is called
    error('Please specify TestPath');
end
path    = contopts.TestPath;
slashes = (path =='/');
slashes = slashes + (path =='\');
sl      = find(slashes);
path    = path(1:sl(end));


%% create data and log files when necessary
datapath = strcat(path,'Data', filesep);
if ~exist(datapath,'dir')
    mkdir(datapath);
    addpath(datapath);
end

logpath  = strcat(path,'Logs', filesep);
if ~exist(logpath,'dir')
    mkdir(logpath);
    addpath(logpath);
end

%% open data file and print header
cds.datapath = datapath;
datafile = strcat(datapath, filename, '.dat');
cds.dataFID  = fopen(datafile,'wt');
fprintf(cds.dataFID, 'CL_MATCONTL2.0\n');
fprintf(cds.dataFID, '%s\n', scriptname);
fprintf(cds.dataFID, '%s\n', prob);
fprintf(cds.dataFID, '%s\n', curve);
fprintf(cds.dataFID, '%f ' , current_time);
fprintf(cds.dataFID, '\n');

%% open log file and print header
logfile = strcat(logpath,'/',filename,'.txt');
if contopts.contL_LogFile
    cds.logFID   = fopen(logfile, 'wt');
else
    cds.logFID = 1;
end

cds.num_sings = 1;
% Header

print_diag(0,'             :  ')
for i = 1:cds.nap
    print_diag(0,'p(%d)              ', cds.ActiveParams(i))
end
if has_period(cds.curve)
%              :  p(2)            norm of point 	 curve function norm 	 step size 
%   1   1   00 :  +5.000000e-01     1.882005e+01     0.000000e+00     1.000000e-02
%   2          :  +5.100000e-01     1.882005e+01     0.000000e+00     1.000000e-02
%   3   2   H  :  +5.128157e-01     1.882005e+01     0.000000e+00     2.815666e-03
  print_diag(0,'period           point norm       curve func norm  step size\n')
else
  print_diag(0,'point norm       curve func norm  step size\n')
end


% load previous bifurcations if we are extending a curve
if exist([fullfile(cds.datapath, cds.runID), '.mat'], 'file')
  if contopts.is_extension
    print_diag(0, 'Loading previously saved bifurcations from file.\n');
    load([fullfile(cds.datapath, cds.runID), '.mat'], 's');
    cds.sout = s;
  else
    print_diag(1, ['A singularities file ' fullfile(cds.datapath, cds.runID) ...
                          '  is is already present, and will be overwritten']);
  end
else
  cds.sout = [];
end

