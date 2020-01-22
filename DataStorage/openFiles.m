function [datafile, logfile] = openFiles()
global cds contopts



prob          = func2str(cds.probfile);
curve         = func2str(cds.curve);

%% create cds.runID
% dbstack returns the function call stack. That is, denote the function that
% called this function by f_1, the function that called f_1 by f_2, the one that
% called f_2 by f_3, and so on. Then dbstack returns an array of structs with
% information about f_i in the i-th place in the array.
stack        = dbstack();
% We expect openFiles to be called by contL, and contL to be called by the
% script that runs the continuation. Hence, we expect stack(3).file to be a
% good name for the files in which the data and logs are saved.
scriptname   = stack(3).file;
% remvome .m from the scriptname
scriptname   = scriptname(1:end-2);
my_timestamp = timestamp();

cds.runID = contopts.Filename;
if isempty(cds.runID) 
  cds.runID = scriptname;
end
% Don't let users accidentally overwrite their data
if contopts.timestamp_in_filenames
  cds.runID = [cds.runID '_' my_timestamp];
end

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
datapath = strcat(path, 'Data', filesep);
if ~exist(datapath,'dir')
    mkdir(datapath);
    addpath(datapath);
end

logpath  = strcat(path, 'Logs', filesep);
if ~exist(logpath,'dir')
    mkdir(logpath);
    addpath(logpath);
end

%% open data file and print header
cds.datapath = datapath;
datafile = strcat(datapath, cds.runID, '.dat');
cds.dataFID  = fopen(datafile,'wt');
fprintf(cds.dataFID, 'CL_MATCONTL2.0\n');
fprintf(cds.dataFID, '%s\n', scriptname);
fprintf(cds.dataFID, '%s\n', prob);
fprintf(cds.dataFID, '%s\n', curve);
fprintf(cds.dataFID, '%s\n', my_timestamp);

%% open log file and print header
logfile = strcat(logpath,'/',cds.runID,'.txt');
if contopts.contL_LogFile
    cds.logFID   = fopen(logfile, 'wt');
else
    cds.logFID = 1;
end

cds.num_sings = 1;


% load previous bifurcations if we are extending a curve
if exist([fullfile(cds.datapath, cds.runID), '.mat'], 'file')
  if contopts.is_extension
    print_diag(0, 'Loading previously saved bifurcations from file.\n');
    load([fullfile(cds.datapath, cds.runID), '.mat'], 's');
    cds.sout = s;
  else
    print_diag(1, ['A singularities file ' fullfile(cds.datapath, cds.runID) ...
                         ' is is already present, and will be overwritten\n']);
  end
else
  cds.sout = [];
end

