function [ failed ] = savePoint( varargin )

% Saves a point to datafile and logfile
% possible arguments:   savePoint(point)
%                       savePoint(point, s)
global cds contopts

point = varargin{1};

if nargin > 2
  error(fprintf(['Wrong number of arguments. Number of arguments' ...
    ' must be 1 or 2. The number of arguments is %d'], nargin));
end
if isequal(cds.curve, @limitcycleL) || ...
   isequal(cds.curve, @single_shooting) || ...
   isequal(cds.curve, @multiple_shooting) 
   print_diag(0,'period: %.6f\n', point.x(end-1));
end
if contopts.every_point_in_separate_mat_file
  % the if statement is to prevent a "Directory already exists" warning
  if ~ exist(fullfile(cds.datapath, cds.runID), 'dir')
    mkdir(cds.datapath, cds.runID)
  end
  % Instead of adding the field savetime to point, we could use modified times
  % of the files from the filesystem, but as far as I know there is no
  % convenient way to read file modified times from matlab.
  point.savetime = clock;
  filename = fullfile(cds.datapath, cds.runID, sprintf('point_%08d', cds.i));
  save(filename, 'point');
end
if nargin == 1
    % 'normal' point
    if cds.i == 1 % first point
        s.index   = 1;
        s.label   = '00';
        s.msg     = 'This is the first point of the curve';
        s.data    = point;
        s.data.ap = cds.ActiveParams;
        s.data.P0 = cds.P0;
        cds.num_sings = 0;
        failed = savePoint(point, s);
        return
    elseif cds.lastpointfound || cds.i >= contopts.MaxNumPoints
        s.index   = cds.i;
        s.label   = '99';
        s.msg     = 'This is last first point of the curve';
        s.data    = point;
        s.data.ap = cds.ActiveParams;
        s.data.P0 = cds.P0;
        
        failed = savePoint(point, s);
        cds.EndTime = clock;
        print_diag(0,'\nElapsed time  = %.1f secs\n', etime(cds.EndTime, cds.StartTime));
        print_diag(0,'Number of Points = %d\n', cds.i);
        
        % save .mat file with singular points and options
        s = cds.sout;
        save([cds.datapath, cds.runID, '.mat'], 's', 'contopts')
        return
    end

    % Carel: print every point on command line
    print_diag(0,'%3d\t%2d\t %s:  '  ,cds.i,[],'  ');
    print_diag(0,'%+.6e  '      ,point.x(end-cds.nap+1:end));
    print_diag(0,'%.6e  '       ,norm(point.x(1:cds.ncoo)));
    print_diag(0,'%.6e  \n'     ,point.R);
    print_diag(0,'Current Step Size: %+e\n',cds.h);
    if (isfield(point,'angle')) % Carel Jonkhout
      % apparently there is no field angle when continuing limit cycles
      print_diag(2,'Angle Between Tangents:  %+e * pi\n',point.angle/pi());
    end
elseif nargin == 2
    % 'singular' point
    s = varargin{2};
    s.index  = cds.i;
    s.data = point;
    s.data.ap = cds.ActiveParams;
    s.data.P0 = cds.P0;
    
    % print to window
    cds.num_sings = cds.num_sings + 1;
    print_diag(0,'%3d\t%2d\t %s:  '  ,cds.i,cds.num_sings,s.label);
    print_diag(0,'%+.6e  '      ,point.x(end-cds.nap+1:end));
    print_diag(0,'%.6e  '       ,norm(point.x(1:cds.ncoo)));
    print_diag(0,'%.6e  \n'     ,point.R);
    try
      if ~ isfield(cds, 'sout'); cds.sout = []; end
      cds.sout = [cds.sout, s];
    catch e
      disp('A minor error occurred in savePoint.m')
      disp(e)
      disp(e.stack)
    end
    if contopts.always_save_s
      s = cds.sout;
      save([cds.datapath, cds.runID, '.mat'], 's', 'contopts')
    end

    
end


if contopts.NewtonPicard
  failed = savePointDataFile(point.x);
elseif isfield(point,'multipliers') %Carel
  failed = savePointDataFile(...
    point.x, point.v, point.h, point.tvals', point.uvals', point.multipliers);
else
  failed = savePointDataFile(...
    point.x, point.v, point.h, point.tvals', point.uvals');
end


function [ failed ] = savePointDataFile( varargin )
% Saves a continuer point to a data file
%   [ failed ] = saveOutput( x,v,h,f,eigs,nsub,asub, etc.... )
% This adds a block to the data file cds.datafile with the folowing
% information
% Line 1: number of arguments passed to savePoint
% Line n: the nth argument to savePoint
global cds
%writes number of lines to save
try
    failed = 0;
    fprintf(cds.dataFID,'%d\n',nargin); %Number data entries to read in a block.
    for i = 1:nargin
        sz = size(varargin{i});
        fprintf(cds.dataFID,'%d %d\n',sz);
        if ~isempty(varargin{i})
          fprintf(cds.dataFID,'%+e ',varargin{i});
          fprintf(cds.dataFID,'\n');
       end
    end
catch
    failed = 1;
end
