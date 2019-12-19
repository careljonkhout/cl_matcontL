function makeEigsPlus()
%   makeEigsPlus:
%       Installs the eigsPlus function for the currenly running release of
%       MATLAB.  Each install is version specific so this command will need
%       to be rerun for each version of MATLAB used.  eigsPlus can also be
%       used on Octave, though not all of its features will be enabled.
%       Running this routine will create up to two .m files in the
%       'private' subfolder of the eigsPlus installation location.
% 
%       MATLAB releases prior to R2012b are NOT supported.
%   
%   NOTE:  
%       To run makeEigsPlus, both the patch and md5 commands must be
%       available on the command line, with the former supporting -u diffs.
%       This should automatically be the case for Linux and Mac users.
%       Windows users may need to install these tools themselves.
%
%   See also cleanEigsPlus, eigsPlus.
%   
%
%   For comments/bug reports, please visit the eigsPlus GitLab webpage:
%   https://gitlab.com/timmitchell/eigsPlus
%
%   eigsPlus Version 2.1, 2016-2018, see AGPL license info below.
%   makeEigsPlus.m introduced in eigsPlus Version 2.0.
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
        error(  'makeEigsPlus:octave',                                  ...
                [   'eigsPlus patching features are only supported on ' ...
                    'MATLAB, not Octave.  eigsPlus may still be used '  ...
                    'on Octave, but not all features will be enabled.'  ]);
    end
    
    % the folder this file is contained in
    base_path       = fileparts(mfilename('fullpath'));
    fs              = filesep();
    diff_path       = [base_path fs getEigsPlusDiffsPath()];
    eigsPlus_name   = [base_path fs 'private'];
      
    mr              = matlabroot;
    sparfun_path    = [mr fs 'toolbox' fs 'matlab' fs 'sparfun'];
    eigs_name       = [sparfun_path fs 'eigs.m'];
    
    % get details for this version of MATLAB (and ensures that it is at
    % least R2012b)
    [suffix,rel_name,newer_matlab] = getEigsPlusSuffix();
    
    % load the saved check_sum variable for the specific version of eigs.m
    load([diff_path fs 'eigs_md5_' suffix '.mat']);

    % get the md5 check sum of eigs, from the running version of MATLAB,
    % and compare to the stored md5 check sum for eigs - they should match!
    [err_status,eigs_cs] = system(['md5 -q ' eigs_name]);
    if err_status
        error(  'makeEigsPlus:checkSumError',                           ...
                [   'Failed to produce a md5 check sum of MATLAB''s '   ...
                    'eigs.m.  Please see the eigsPlus webpage for '     ...
                    'support options.'  ]                               );
    end
    if any(check_sum ~= eigs_cs)
        if newer_matlab
            writeConfigurationFile(rel_name,true);
            error('makeEigsPlus:unsupportedMatlab',                     ...
                    [   'Newer unsupported MATLAB version detected. '   ...
                        'eigsPlus may still be used, but not all '      ...
                        'features will be enabled.  Please check the '  ...
                        'eigsPlus website if a newer version of '       ...
                        'eigsPlus is available.'    ]                   );      
        end
        error(  'makeEigsPlus:checkSumMismatch',                        ...
                [   'The md5 check sum of MATLAB''s eigs.m does not '   ...
                    'match the known version.  Please see the eigsPlus '...
                    'webpage for support options.'  ]                   );
    end
     
    % the filename for the particular patch file
    patch_name      = [diff_path fs 'eigs_patch_' suffix '.txt'];
    
    % the filename for the eigsPlus patched version of eigs to be created
    eigsPlus_name   = [eigsPlus_name fs getEigsPlusName(rel_name)];
    
    % the commandline string 
    patch_cmd       = sprintf('patch -s -u %s -i %s -o %s', eigs_name,  ...
                                                            patch_name, ...
                                                            eigsPlus_name);
                                                        
    % attempt to apply the patch
    [err_status] = system(patch_cmd);                                   
    if err_status
        error(  'makeEigsPlus:patchError',                              ...
                [   'Failed to install eigsPlus patch.  Please see '    ...
                    'the eigsPlus webpage for support options.' ]       );
    end
    
    writeConfigurationFile(rel_name);
    
    fprintf('eigsPlus successfully installed for MATLAB R%s.\n',rel_name);
end