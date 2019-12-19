function eigsPlus_name = getEigsPlusName(rel_name)
%   getEigsPlusName:
%       Returns the version-specific filename of the patched eigsPlus
%       routine.
%
%   INPUT:
%       rel_name        'XXXXX' string for the Matlab release (e.g. 2016b)
%
%   OUTPUT:
%       eigsPlus_name   'eigsPlus_RXXXXX.m' string (e.g. eigsPlus_R2016b.m)  
%
%
%   For comments/bug reports, please visit the eigsPlus GitLab webpage:
%   https://gitlab.com/timmitchell/eigsPlus
%
%   getEigsPlusName.m introduced in eigsPlus Version 2.0.
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

    eigsPlus_name = ['eigsPlus_R' rel_name '.m'];
end