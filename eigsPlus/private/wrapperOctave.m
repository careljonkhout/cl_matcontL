function varargout = wrapperOctave(varargin)
%   wrapperEigsPlus:
%       This wrapper routine modifies the output arguments of Octave's
%       version of the eigs routine.
%
%       First, it always returns the eigenvalues as a column vector d of
%       entries, as opposed to Matlab's eigs, which uses a column vector d
%       when only the eigenvalues are requested and then a diagonal matrix
%       D when the eigenvectors V are additionally requested.
% 
%       Second, it removes any unconverged Ritz values (and the
%       corresponding eigenvector approximations) from the eigenvalue and
%       eigenvector output arguments d and V.
%
%       eigsPlus makes use of wrapper functions for compatibility reasons,
%       so that the eigsPlus output argument format can also be used with
%       Octave as well as Matlab's default eigs routine.  It also means
%       that less modifications to eigs are needed, which simplifies
%       maintaining eigsPlus.
%
%   INPUT:   
%       varargin        normal input arguments to eigs
%   
%   OUTPUT:
%       d               = wrapperOctave(...)
%       [V,d]           = wrapperOctave(...)
%       [V,d,flag]      = wrapperOctave(...)
%       
%       d               Column vector of p converged eigenvalues.  Can be
%                       empty if no eigenvalues converged.  eigsPlus always
%                       returns d as a column vector, regardless of whether
%                       or not the eigenvectors are also requested.
% 
%       V               Matrix of p eigenvectors corresponding to the p
%                       converged eigenvalues given in d.  Can be empty if
%                       no eigenvalues converged.
% 
%       flag            Default 'flag' output argument of eigs.
% 
%
%   For comments/bug reports, please visit the eigsPlus GitLab webpage:
%   https://gitlab.com/timmitchell/eigsPlus
%
%   wrapperOctave.m introduced in eigsPlus Version 2.0.
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

    % We need to get the eigenvectors no matter what, since for Octave, we
    % need them to help determine which eigenvalues did NOT converge.
    [V,D,flag]      = eigs(varargin{:});
    d               = diag(D);
    
    % Octave uses zeros for unconverged Ritz values, for the eigenvalues
    % and the corresponding eigenvectors.  Thus, to remove unconverged Ritz
    % values from Octave, we can check both that both returned eigenvalues
    % and eigenvectors are zero.  We cannot just assume that the zero
    % eigenvalues actually represent unconverged Ritz values since zero
    % itself could be an actual eigenvalue of the matrix.
    
    % Note that in Matlab R2013b and later, eigs uses NaNs as placeholders
    % for unconverged Ritz values, but it only does so for the eigenvalues.
    % If eigenvectors are also requested, Matlab's eigs will still return
    % approximations to the eigenvectors for the these unconverged Ritz
    % values.  Either way, in Matlab, it is sufficient to just check the
    % eigenvalues for NaNs, which is not only simpler but also quicker.
    
    % I have provided Octave with a bug report about this Matlab
    % compatility issue.  This issue exists in at least Octave 4.0.0.  
    % It will hopefully be addressed in the 4.4.0 release. 

    % If not all the requested eigenvalues converged, we need to eliminate
    % the unconverged values.
    if flag
        % get the indices for the converged eigenvalues
        indx        = any(V,1) | any(d,2)'; 
        % replace the unconverged entries with NaNs
        d(~indx)    = nan;
    end
   
    if nargout < 2
        [varargout{1:nargout}] = d;
    else
        out_args = {V,d,flag};
        [varargout{1:nargout}] = out_args{1:nargout};
    end
end