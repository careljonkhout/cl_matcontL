function makeEigsPlusDiffs()
%   makeEigsPlusDiffs:
%       This routine is intended for eigsPlus developer use only.
% 
%       It makes all the version-specific eigsPlus patches and md5 check
%       sums.
%
%       eigsPlus only ships diff patches to the various versions of
%       Matlab's eigs routine, as including the complete, eigsPlus-modified
%       versions would violate The MathWorks copyright on their source
%       code.
%   
%       To run this routine, one must have have valid access to the
%       original eigs source code and to a modified version which includes
%       the eigsPlus modifications.  A version-specific diff patch will be
%       created for each pair of the original source and its corresponding
%       eigsPlus-modified source file.
%
%       NOTE: getEigsPlusSuffix.m must also be updated when adding support 
%       for new releases of MATLAB.
%   
%   REQUIREMENTS:
%   
%   1)  In the parent directory of the 'private' folder containing this
%       file, there must exist a 'sources' folder (i.e. 'sources' and
%       'private' exist in the same folder, namely the eigsPlus
%       installation folder).  The 'sources' folder must contain two
%       subfolders: 'matlab' and 'modified'.
%
%       The 'matlab' folder must contain subfolders 'RXXXXX' (e.g. R2016b)
%       for every release of Matlab supported by eigsPlus.  Inside each of
%       these subfolders must be the original Matlab eigs.m source file for
%       that particular version of Matlab.
%
%       The 'modified" folder must contain each eigsPlus modified version
%       of the version-specific Matlab eigs routine.  Since eigs is not
%       always modified with every new release of Matlab, each eigsPlus
%       version may pertain to a range of Matlab releases.  For examples:
% 
%           eigsPlus_R2015a_R2016a.m
%           eigsPlus_R2016b.m
%   
%   2)  The diff and md5 commands must be available on the system terminal,
%       with the former supporting the -u option.  These commmands should
%       be available by default on all Linux and Mac operating systems.
%       Windows users may need to install them.
%      
%
%   For comments/bug reports, please visit the eigsPlus GitLab webpage:
%   https://gitlab.com/timmitchell/eigsPlus
%
%   makeEigsPlusDiffs.m introduced in eigsPlus Version 2.0.
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
    
    diff_type   = 'eigsPlus_';
    pre_chars   = length(diff_type);

    % this is the parent folder (since this is the private folder) 
    base_path   = fileparts(fileparts(mfilename('fullpath')));
    fs          = filesep();
    org_path    = [base_path fs 'sources' fs 'matlab'];
    mods_path   = [base_path fs 'sources' fs 'modified'];
    diff_path   = [base_path fs getEigsPlusDiffsPath()];
    
    if ~exist(org_path,'dir') || ~exist(mods_path,'dir')
        error(  'makeEigsPlusDiffs:sourceFilesMissing',                 ...
                [   'Source files for eigs.m and the corresponding '    ...
                    'eigsPlus versions are missing.  Cannot make '      ...
                    'eigsPlus patches and md5 check sums.'  ]           );
    end
    
    if ~exist(diff_path,'dir')
        mkdir(diff_path);
    end
    
    eigsPlus_mods       = dir([mods_path fs diff_type '*.m']);
    
    for j = 1:length(eigsPlus_mods)
        eigsPlus_mod    = eigsPlus_mods(j).name;
        [~,name,ext]    = fileparts(eigsPlus_mod);  
        
        % skip files that are .m or being with the diff_type prefix
        if ~strcmpi('.m',ext)
            continue
        end
        if ~strncmpi(name,diff_type,pre_chars)
            continue
        end
        
        patch_rel       = name(1+pre_chars:end);
        patch_rels      = strsplit(patch_rel,'_');
        
        org_name        = [org_path fs patch_rels{1} fs 'eigs.m'];
        patched_name    = [mods_path fs eigsPlus_mod];
        diff_name       = [diff_path fs 'eigs_patch_' patch_rel '.txt'];
        md5_name        = [diff_path fs 'eigs_md5_' patch_rel '.mat'];
        
        [~,check_sum]   = system(['md5 -q ' org_name]);
        save(md5_name,'check_sum');
        
        % need to diff to our stored originals, not my version of matlab
        system(['diff -u ' org_name ' ' patched_name ' > ' diff_name]);
    end
    
end