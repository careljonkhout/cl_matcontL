function [ varargout ] = loadfile( filename, varargin )
% Loads a continuer point from a data file
%   [ failed ] = saveOutput( x,v,h,f,eigs,nsub,asub, etc.... )
% This adds a block to the data file cds.datafile with the folowing
% information
% Line 1: number of arguments passed to savePoint
% Line n: the nth argument to savePoint
%writes number of lines to save
FID = fopen(filename,'r');
for i=1:5
    fgetl(FID);
end
% fprintf(1,'%s\n',fgetl(FID));
% fprintf(1,'Testrun: %s\n',fgetl(FID));
% fprintf(1,'Problem: %s\n',fgetl(FID));
% fprintf(1,'Curve:   %s\n',fgetl(FID));
% date = cell2mat(textscan(fgetl(FID),'%f'));
% 
% fprintf(1,'Date:    %d-%d-%d %02d:%02d:%02.f\n',date);


j=1;
while 1
    tline = fgetl(FID);
    if ~ischar(tline), break, end
    blocksize = cell2mat(textscan(tline,'%f'));
    for i=1:blocksize
        tline = fgetl(FID);
    
        arraydim = cell2mat(textscan(tline,'%f'));
        nextread = zeros(arraydim');
        for k = 1:arraydim(2)
            tline = fgetl(FID);
            temp = cell2mat(textscan(tline,'%f'));
            nextread(:,k) = temp;
        end
        
        if arraydim(2) == 1
            varargout{i}(1:arraydim(1),j) = nextread;
        elseif arraydim(2) > 1
            varargout{i}{j} = nextread';
        end
    end
    j=j+1;
end
fclose(FID);