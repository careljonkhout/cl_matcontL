function varargout = wrapperEigsPlus(eigsPlus_fn,varargin)
%   wrapperEigsPlus:
%       This wrapper routine modifies the output arguments of eigsPlus (all
%       versions) and/or Matlab's default eigs (only R2013b or later are
%       supported).
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
%       eigsPlus makes use of a wrapper functions for compatibility
%       reasons, so that the eigsPlus output argument format can also be
%       used with Octave as well as Matlab's default eigs routine.  It also
%       means that less modifications to eigs are needed, which simplifies
%       maintaining eigsPlus.
%
%   INPUT:
%       eigsPlus_fn     function handle to an eigsPlus or default eigs 
%                       routine
%   
%       varargin        normal input arguments to eigs
% 
%   OUTPUT:
%       d               = wrapperEigsPlus(...)
%       [V,d]           = wrapperEigsPlus(...)
%       [V,d,iters]     = wrapperEigsPlus(...)
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
%       iters           The number of ARPACK iterations incurred (eigsPlus
%                       only; for Matlab's default eigs, this is the 'flag'
%                       output argument).
%     
%
%   For comments/bug reports, please visit the eigsPlus GitLab webpage:
%   https://gitlab.com/timmitchell/eigsPlus
%
%   wrapperEigsPlus.m introduced in eigsPlus Version 2.0.
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


    % Process the output arguments to the eigsPlus standard

    % Calls the version specific eigsPlus_R20XXX.m
    [varargout{1:nargout}] = eigsPlus_fn(varargin{:});

    % Note that in Matlab R2013b and later, eigs uses NaNs as placeholders
    % for unconverged Ritz values, but it only does so for the eigenvalues.
    % If eigenvectors are also requested, Matlab's eigs will still return
    % approximations to the eigenvectors for the these unconverged Ritz
    % values.  Either way, in Matlab, it is sufficient to just check the
    % eigenvalues for NaNs, which is not only simpler but also quicker.

    if nargout < 2
        d               = varargout{1};
        indx            = ~isnan(d);
        varargout{1}    = d(indx);
    else
        V               = varargout{1};
        d               = diag(varargout{2});
        indx            = ~isnan(d);
        varargout{1}    = V(:,indx);
        varargout{2}    = d(indx);
    end
end