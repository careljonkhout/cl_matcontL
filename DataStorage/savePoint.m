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

  % the if statement is to prevent a "Directory already exists" warning
  if ~ exist(fullfile(cds.datapath, cds.runID), 'dir')
    mkdir(cds.datapath, cds.runID)
  end
  point.savetime = clock;
  filename = fullfile(cds.datapath, cds.runID, sprintf('point_%08d', cds.i));
  save(filename, 'point');

  if nargin == 1
      % 'normal' point
      if cds.first_point % first point
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
          print_diag(0,'\nElapsed time  = %.1f secs\n', ...
                        etime(cds.EndTime, cds.StartTime));
          print_diag(0,'Number of Points = %d\n', cds.i);

          % save .mat file with singular points and options
          s = cds.sout;
          save([cds.datapath, cds.runID, '.mat'], 's', 'contopts')
          return
      end

      label = '  ';
      is_special = false;
      print_point(point, label, is_special);
      if isfield(point, 'angle')
        print_diag(2, 'Angle Between Tangents:  %+f degrees\n', ...
                point.angle / pi * 180);
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
      is_special = true;
      print_point(point,s.label,is_special);
      try
        if ~ isfield(cds, 'sout'); cds.sout = []; end
        cds.sout = [cds.sout, s];
      catch e
        disp('A minor error occurred in savePoint.m')
        disp(e)
        disp(e.stack)
      end
      
      if ~ isempty(contopts.singularity_callback)
        contopts.singularity_callback(s)
      end

      s = cds.sout;
      save([cds.datapath, cds.runID, '.mat'], 's', 'contopts')
  end


  if isfield(point,'multipliers')
    failed = savePointDataFile(...
      point.x, point.v, point.h, point.tvals', point.uvals', point.multipliers);
  else
    failed = savePointDataFile(...
      point.x, point.v, point.h, point.tvals', point.uvals');
  end
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
end




function print_point(point, label, is_special)
  global cds contopts

  if ~ isfield(cds, 'n_points_printed')
    cds.n_points_printed = 0;
  end
 
  if mod(cds.n_points_printed, 20) == 0 || contopts.console_output_level >= 1
    PrintOutputCaption()
  end
  
  if is_special
    print_diag(0, '%3d  %2d  %4s :  ' , cds.i ,cds.num_sings, label);
    priority = 0;
  else
    print_diag(1, '%3d           :  '  , cds.i                     );
    priority = 1;
  end
 
% example of output produced by this function:
%  S  ID   PT:     p(2)          norm of point   curve function norm 	 step size 
%  1	 1   00:  +5.000000e-01    6.760362e+01    0.000000e+00          1.000000e-02    
%  2         :  +5.100000e-01    6.760362e+01    0.000000e+00          1.000000e-02  

  print_diag(priority, '%+.6e     ', point.x(end-cds.nap+1:end));
  if has_period(cds.curve)
    print_diag(priority, '%.6e     ', point.x(end-1));
  end
  print_diag(priority, '%.6e     ', norm(point.x(1:cds.ncoo)));
  print_diag(priority, '%.6e     ', point.R);
  print_diag(priority, '%12.6e'   , point.h);
  
  if has_period(cds.curve) && has_nonempty_field(point, 'multipliers')
    distance_to_one     = abs(point.multipliers - 1);
    log_10_of_deviation = log10(min(distance_to_one));
    print_diag(priority, '%10.3f', -log_10_of_deviation);
    m = point.multipliers(1);
    if abs(imag(m)) > cds.deviation_of_trivial_multiplier
      print_diag(priority, '%18.6f +/- %.6fi norm: %.6f', ...
                        real(m), abs(imag(m)), abs(m));
    else
      print_diag(priority, '%18.6f', real(m));
    end
  end
  print_diag(priority, '\n');
  if is_special
    cds.n_points_printed = cds.n_points_printed + 1;
  end
end

function has_nonempty_field = has_nonempty_field(struct, field)
  has_nonempty_field = isfield(struct, field) && ~ isempty(struct.(field));
end

  