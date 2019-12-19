function cleanEigsPlusAll()
%   cleanEigsPlusAll:
%       This routine removes all the version-specific files of eigsPlus
%       that are currently installed.  These are the files that
%       makeEigsPlus.m creates in this 'private' folder when it is run.
%
%
%   For comments/bug reports, please visit the eigsPlus GitLab webpage:
%   https://gitlab.com/timmitchell/eigsPlus
%
%   cleanEigsPlusAll.m introduced in eigsPlus Version 2.0.
%
% =========================================================================
% |  eigsPlus                                                             |
% |  Copyright (C) 2016 Tim Mitchell                                      |
% |                                                                       |
% |  This file is part of eigsPlus.                                       |
% |                                                                       |
% |  eigsPlus is free software: you can redistribute it and/or modify     |
% |  it under the terms of the GNU Affero General Public License as       |
% |  published by the Free Software Foundation, either version 3 of       |
% |  the License, or (at your option) any later version.                  |
% |                                                                       |
% |  eigsPlus is distributed in the hope that it will be useful,          |
% |  but WITHOUT ANY WARRANTY; without even the implied warranty of       |
% |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        |
% |  GNU Affero General Public License for more details.                  |
% |                                                                       |
% |  You should have received a copy of the GNU Affero General Public     |
% |  License along with this program.  If not, see                        |
% |  <http://www.gnu.org/licenses/>.                                      |
% =========================================================================

    base_path           = fileparts(mfilename('fullpath'));
    potential_files     = dir([base_path filesep '*_R*.m']);
    n_files             = length(potential_files);
    file_pattern        = '(eigsPlus|useDefaultEigs)_R\d{4}[ab]{1}\.m';
    version_pattern     = '\d{4}[ab]{1}';
    versions            = cell(1,n_files);
    count               = 0;
        
    % find the versions that are installed
    for j = 1:n_files     
        name    = potential_files(j).name;
        n_chars = length(name);
        
        if n_chars ~= 17 && n_chars ~= 23
            continue
        end
        
        indx = regexp(name,file_pattern);
        if isempty(indx) || indx(1) > 1
            continue
        end
        
        % we are now assured that it is an eigsPlus or config file
        % extract the version number from the file name
        indx            = regexp(name,version_pattern);
        % add it to the list
        count           = count + 1;
        versions{count} = name(indx:indx+4);
    end
    
    % crop empty entries and get rid of the duplicates
    versions    = unique(versions(1:count));
    
    for j = 1:length(versions)
        cleanEigsPlusVersion(versions{j});
    end
end