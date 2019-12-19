function cleanEigsPlus(varargin)
%   cleanEigsPlus:
%       Removes the version-specific eigsPlus files that are created
%       anytime makeEigsPlus is run.
%
%   USAGE:
%       % Removes all versions of the files created by each call to
%       % makeEigsPlus(), called on whichver versions of MATLAB.
%       > cleanEigsPlus()
%
%       % Remove  a specific version of the files created by makeEigsPlus()
%       % The single input must be a string with the following form:
%       %
%       %       XXXX[ab] where X is a digit 0-9, e.g. 2015a
%       %   
%       % and the string/release is at least 2006a
%       > cleanEigsPlus('2015a')
%
%   NOTE:  
%       This command is only applicable on MATLAB, not Octave, as the 
%       eigs patching features are MATLAB-only features.
%
%   See also eigsPlus, makeEigsPlus.
%   
%
%   For comments/bug reports, please visit the eigsPlus GitLab webpage:
%   https://gitlab.com/timmitchell/eigsPlus
%
%   eigsPlus Version 2.1, 2016-2018, see AGPL license info below.
%   cleanEigsPlus.m introduced in eigsPlus Version 2.0.
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

    if isOctave()
        error(  'cleanEigsPlus:octave',                                 ...
                [   'Nothing to clean.\n'                               ...
                    'eigsPlus patching features are only supported on ' ...
                    'MATLAB, not Octave.  eigsPlus may still be used '  ...
                    'on Octave, but not all features will be enabled.'  ]);
    end
    
    if nargin < 1
        cleanEigsPlusAll();
    else
        version_str = varargin{1};
    
        assertMatlabReleaseStr(version_str);
        version_str = lower(version_str);
        
        % make sure user's version is at least 2012b, the oldest version 
        % supported by eigsPlus
        if strCompareTo('2012b',version_str) > 0
            error(  'cleanEigsPlus:unsupportedMatlab',                  ...
                    'eigsPlus requires MATLAB R2012b or newer.'         );             
        end
        
        cleanEigsPlusVersion(version_str);
    end
end