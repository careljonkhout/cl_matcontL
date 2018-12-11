function [datafile, logfile] = openFiles()
global cds contopts

%% create filename
[ST,I]   = dbstack();
testname = ST(end).file;
testname = testname(1:end-2);
prob     = func2str(cds.probfile);
curve    = func2str(cds.curve);
thetime  = clock;
strtime  = sprintf('%d_%02d_%02d',thetime(4:end-1),round(thetime(end)));

if isempty(contopts.Filename)
    fname    = strcat(testname,'_',date(),'_',strtime);
else
    fname = contopts.Filename;
end
cds.runID = fname;

%% Find path
if isempty(contopts.TestPath)   % note that it is attempted to set testpath when contset() is called
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
datafile = strcat(datapath,fname,'.dat');
cds.dataFID  = fopen(datafile,'wt');
fprintf(cds.dataFID,'MATCONTL2.0\n');
fprintf(cds.dataFID,'%s\n',testname);
fprintf(cds.dataFID,'%s\n',prob);
fprintf(cds.dataFID,'%s\n',curve);
fprintf(cds.dataFID,'%f ',thetime);
fprintf(cds.dataFID,'\n');

%% open log file and print header
logfile = strcat(logpath,'/',fname,'.txt');
if contopts.contL_LogFile
    cds.logFID   = fopen(logfile, 'wt');
else
    cds.logFID = 1;
end

cds.num_sings = 1;
% Header
print_diag(0,' S\tID\t PT:  ')
for i = 1:cds.nap
    print_diag(0,'p(%d)            ', cds.ActiveParams(i))
end
print_diag(0,'||u||          ||f||\n')