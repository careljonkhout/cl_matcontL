function cleanEigsPlusVersion(rel_name)
%   cleanEigsPlusVersion:
%       This routine deletes the version-specific files of eigsPlus that
%       are currently installed for the version of Matlab specified by
%       rel_name. These are the files that makeEigsPlus.m creates in this
%       'private' folder when it is run on the corresponding version of
%       Matlab.
% 
%   INPUT:
%       rel_name    string specifying the Matlab release as XXXX[ab]
%                   (e.g. '2016b')
%    
%
%   For comments/bug reports, please visit the eigsPlus GitLab webpage:
%   https://gitlab.com/timmitchell/eigsPlus
%
%   cleanEigsPlusVersion.m introduced in eigsPlus Version 2.0.
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

    base_path       = fileparts(mfilename('fullpath'));
    fs              = filesep();
    eigsPlus_file   = [base_path fs getEigsPlusName(rel_name)];
    config_file     = [base_path fs getUseDefaultEigsName(rel_name)];
    
    if exist(eigsPlus_file,'file')
        delete(eigsPlus_file);
    end
    
    if exist(config_file,'file')
        delete(config_file);
    end
end