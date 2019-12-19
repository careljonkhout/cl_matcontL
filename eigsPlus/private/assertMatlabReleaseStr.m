function assertMatlabReleaseStr(s)
%   assertMatlabReleaseStr:
%       Throws an assertion error if s is not a string of 5 chars in the 
%       following format used to denote Matlab release names:
%
%           XXXX[ab] where X is a digit 0-9, e.g. 2015a
%          
%       and the string/release is at least 2006a
%       
%   USASE:
%       assertMatlabReleaseStr('2015a')     % succeeds
%       assertMatlabReleaseStr('R2015a')    % throws an error
%
%
%   For comments/bug reports, please visit the eigsPlus GitLab webpage:
%   https://gitlab.com/timmitchell/eigsPlus
%
%   assertMatlabReleaseStr.m introduced in eigsPlus Version 2.0.
%
% =========================================================================
% |  assertMatlabReleaseStr.m                                             |
% |  Copyright (C) 2016 Tim Mitchell                                      |
% |                                                                       |
% |  This file is originally from URTM.                                   |
% |                                                                       |
% |  URTM is free software: you can redistribute it and/or modify         |
% |  it under the terms of the GNU Affero General Public License as       |
% |  published by the Free Software Foundation, either version 3 of       |
% |  the License, or (at your option) any later version.                  |
% |                                                                       |
% |  URTM is distributed in the hope that it will be useful,              |
% |  but WITHOUT ANY WARRANTY; without even the implied warranty of       |
% |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        |
% |  GNU Affero General Public License for more details.                  |
% |                                                                       |
% |  You should have received a copy of the GNU Affero General Public     |
% |  License along with this program.  If not, see                        |
% |  <http://www.gnu.org/licenses/>.                                      |
% =========================================================================
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

    msg = 'Version must be specified in XXXX[ab] format, e.g. ''2015a''';
    assert(ischar(s) && length(s) == 5, msg);
    
    % let's not care whether or not uppercase or lower case was used for 
    % A or B.
    s   = lower(s);
    y   = str2double(s(1:4));
    ab  = s(end);
    assert(ab == 'a' || ab == 'b', msg);
    assert(~isnan(y) && isAnInteger(y) && y >= 2006, msg);
end