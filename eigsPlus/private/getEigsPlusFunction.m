function fn = getEigsPlusFunction()
%   getEigsPlusFunction:
%       Returns a function handle to the version-specific patched eigsPlus
%       routine.  This will throw an error if eigsPlus is not installed for
%       the currently running version of MATLAB or if eigsPlus is
%       improperly configured (the latter which only applies to 2013a or
%       earlier releases).
%
%
%   For comments/bug reports, please visit the eigsPlus GitLab webpage:
%   https://gitlab.com/timmitchell/eigsPlus
%
%   getEigsPlusFunction.m introduced in eigsPlus Version 2.0.
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
    
    % the folder this file is contained in
    base_path       = fileparts(mfilename('fullpath'));
    
    [~,rel_name]    = getEigsPlusSuffix();
    
    % filenames, with .m suffix
    config_fn       = getUseDefaultEigsName(rel_name);
    eigsPlus_fn     = getEigsPlusName(rel_name); 
    
    % fullpath to files
    fs              = filesep();
    config_name     = [base_path fs config_fn];
    eigsPlus_name   = [base_path fs eigsPlus_fn];
    
    % make sure eigsPlus is installed for this version of MATLAB
    if ~exist(config_name,'file') || ~exist(eigsPlus_name,'file')
        error(  'eigsPlus:notInstalled',                            ...
                [   'New version of MATLAB detected, run '          ...
                    'makeEigsPlus() to install/update eigsPlus.']   );
    end
    
    % to make a function handle, we need to lop off the .m on string
    use_eigs_fn     = str2func(config_fn(1:end-2));
    
    if use_eigs_fn()
        if ~atLeastVersion('2013b')
            error(  'eigsPlus:invalidConfiguration',                    ...
                    [   'eigsPlus can only be downgraded to MATLAB''s ' ...
                        'default version of eigs on R2013b and later.'  ]);
        end
        eigs_fn     = @eigs;
    else
        % again, lop off the .m first! 
        eigs_fn     = str2func(eigsPlus_fn(1:end-2));
    end
    fn              = @(varargin) wrapperEigsPlus(eigs_fn,varargin{:});
end