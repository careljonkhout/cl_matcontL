function [ failed ] = savePoint( varargin )

% Saves a point to datafile and logfile
% possible arguments:   savePoint(point)
%                       savePoint(point, s)
global cds contopts

point = varargin{1};
if nargin == 1
    % 'normal' point
    if cds.i == 1 % first point
        s.index   = 1;
        s.label   = '00';
        s.msg     = 'This is the first point of the curve';
        s.data    = point;
        s.data.ap = cds.ActiveParams;
        s.data.P0 = cds.P0;
        
        cds.sout = [];
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
    
    print_diag(1,'%3d\t%2d\t %s:  '  ,cds.i,[],'  ');
    print_diag(1,'%+.6e  '      ,point.x(cds.ncoo+1:cds.ncoo+cds.nap));
    print_diag(1,'%.6e  '       ,norm(point.x(1:cds.ncoo)));
    print_diag(1,'%.6e  \n'     ,point.R);
    print_diag(1,'Current Step Size: %+e\n',cds.h);
    if (isfield(point,'angle')) % Carel Jonkhout
      % apparently there is no field angle when continuing limit cycles
      print_diag(2,'Angle Between Tangents:  %+e * pi\n',point.angle/pi());
    end
    failed = savePointDataFile(point.x, point.v, point.h, point.tvals', point.uvals');
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
    print_diag(0,'%+.6e  '      ,point.x(cds.ncoo+1:cds.ncoo+cds.nap));
    print_diag(0,'%.6e  '       ,norm(point.x(1:cds.ncoo)));
    print_diag(0,'%.6e  \n'     ,point.R);
    failed = savePointDataFile(point.x, point.v, point.h, point.tvals', point.uvals');
    cds.sout = [cds.sout, s];
    
%     if cds.i == contopts.Cont_MaxNumPoints && ~cds.lastpointfound, 
%         cds.lastpointfound = 1; 
%         savePoint(point); 
%     end
else
    error('Wrong number of arguments')
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
    fprintf(cds.dataFID,'%d\n',nargin); %Number data entried to read in a block.
    for i = 1:nargin
        sz = size(varargin{i});
        fprintf(cds.dataFID,'%d %d\n',sz);
        for j= 1:sz(2)
            fprintf(cds.dataFID,'%+e ',varargin{i});
            fprintf(cds.dataFID,'\n');
        end
    end
catch
    failed = 1;
end
