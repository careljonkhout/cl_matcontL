function [eigsPlus_suffix,rel_name,newer_matlab] = getEigsPlusSuffix()
%   getEigsPlusSuffix:
%       Gets the necessary eigsPlus and Matlab version data for the
%       currently running version of Matlab.
% 
%   OUTPUT:
%       eigsPlus_suffix     'RXXXXX' or 'RXXXXX_RXXXXX' string used as the 
%                           suffix for each eigsPlus patch and md5 check
%                           sum of the Matlab's eigs routine.  These
%                           patches and check sums either pertain to a
%                           single version of Matlab or a range of them.
%
%                           Examples:
%                               eigs_patch_R2015a_R2016a.txt
%                               eigs_patch_R2016b.txt
%                               eigs_md5_R2016b.mat
%
%       rel_name            'XXXXX' string for the Matlab release 
%                           (e.g. 2016b)
%
%       newer_matlab        logical indicating whether the currently 
%                           running release of Matlab is newer than the
%                           latest explicitly-supported eigsPlus release.
% 
%
%   For comments/bug reports, please visit the eigsPlus GitLab webpage:
%   https://gitlab.com/timmitchell/eigsPlus
%
%   getEigsPlusSuffix.m introduced in eigsPlus Version 2.0.
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

    newer_matlab    = false;
    rel_name        = lower(version('-release'));
    
    switch rel_name
        case '2012b'
            eigsPlus_suffix = 'R2012b';
        case '2013a'
            eigsPlus_suffix = 'R2013a';
        case '2013b'
            eigsPlus_suffix = 'R2013b';
        case {'2014a','2014b'}
            eigsPlus_suffix = 'R2014a_R2014b';
        case {'2015a','2015b','2016a'}
            eigsPlus_suffix = 'R2015a_R2016a';
        case {'2016b'}
            eigsPlus_suffix = 'R2016b';
        case {'2017a'}
            eigsPlus_suffix = 'R2017a';
        case {'2017b'}
            eigsPlus_suffix = 'R2017b';
        case {'2018a'}
            eigsPlus_suffix = 'R2018a';
        otherwise
            if ~atLeastVersion('2012b')
                error(  'eigsPlus:unsupportedMatlab',                   ...
                        'eigsPlus requires Matlab R2012b or newer.'     );
            else
                newer_matlab = true;
            end
            eigsPlus_suffix = 'R2018a';
    end
end