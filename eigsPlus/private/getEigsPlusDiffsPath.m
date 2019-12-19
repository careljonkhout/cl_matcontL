function eigsPlus_diffs_path = getEigsPlusDiffsPath()
%   getEigsPlusDiffsPath:
%       Returns the name of the folder to contain the version-specific diff
%       patches for creating a patch versions of Matlab's eigs routine.
%
%
%   For comments/bug reports, please visit the eigsPlus GitLab webpage:
%   https://gitlab.com/timmitchell/eigsPlus
%
%   getEigsPlusDiffsPath.m introduced in eigsPlus Version 2.0.
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

    eigsPlus_diffs_path  = 'diffs';
end