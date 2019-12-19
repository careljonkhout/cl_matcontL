function varargout = eigsPlus(varargin)
%   eigsPlus:
%       eigsPlus is an enhanced version of MATLAB's eigs command that
%       provides the following modifications:
%
%       1)  Only the converged eigenvalues and eigenvalue eigenvector pairs 
%           are returned, so the user no longer needs to worry about
%           accidentally using the placeholder zeros (R2013a and earlier)
%           or NaNs for unconverged Ritz values in a subsequent
%           computation.  If none converge, empty arrays are returned.
%
%       2)  The eigenvalues are always returned in a column vector, and 
%           never a diagonal matrix, unlike MATLAB's eigs routine which
%           will return the eigenvalues as a diagonal matrix when
%           eigenvectors are always requested.
%
%       3)  If no eigenvalues are resolved, on R2016b and earlier, eigs
%           would throw an error; on these versions, eigsPlus will instead
%           throw a warning. On R2017a, eigs also adopted the same
%           convention of throwing a warning instead of an error. However,
%           on R2017b, eigs no longer throws a warning or an error. To
%           maintain consistency, as of R2017b, eigsPlus also no longer
%           throws a warning. Recall in the case that no Ritz value
%           converges, the returned eigenvalues and eigenvectors will be
%           empty.
% 
%       4)  The third output argument, flag, has been repurposed to return
%           the number of iterations that is incurred; the user can always
%           equivalently check if the number of eigenvalues returned is
%           less than what the user requested by checking the dimensions of
%           the eigenvalue and/or eigenvector output arguments. For R2017a
%           and earlier, this is the number of ARPACK iterations. On newer
%           releases, this is the number of Krylov Schur iterations but
%           note that for problems with symmetric positive definite B
%           matrices using shift-invert, the Krylov Schur method may be
%           restarted once, if the first fails; in this case, the number of
%           iterations is the sum of two attempts.
%
%       5)  R2017a and earlier only: the user may set opts.isreal to false
%           to force the complex ARPACK routines to be used even if matrix
%           A is real and opts.v0 is either real or complex (whereas
%           MATLABs's eigs would throw an error in these cases).  eigsPlus
%           will throw an error in either of these cases if opts.isreal is
%           either not provided or if it is set to true.
%
%   NOTE:   On Octave, only features (1) and (2) are enabled. 
%
%   NOTE:   eigsPlus may be slightly slower the first time it is called.
%           If you are running sensitive timing experiments, ensure that 
%           eigsPlus has been called at least once before commencing timing 
%           evaluations.  
%
%   eigsPlus was first created to handle a problem that can arise with eigs
%   on MATLAB (R2013a and earlier) and on current versions of Octave (4.0).
%   If a user requests k eigenvalues but ARPACK then fails to resolve k
%   eigenvalues, then, on the affected platforms, eigs will return the
%   converged eigenvalues along with zeros as placeholders for the
%   remaining unconverged Ritz values.  This creates at least two problems.
%
%   First, it can be difficult to remove the placeholder zeros from the
%   list of converged eigenvalues, particularly if zero is actually an
%   eigenvalue of the matrix in question.
%
%   Second, the user may unintentionally use these placeholder zeros in
%   subsequent computations.  For example, if the user is interested in
%   computing the spectral abscissa of a stable matrix A (continuous-time),
%   they may write the following code:
%   
%       max(real(eigs(A,6,'LR')))
%
%   However, if less than 6 Ritz values converge, the result will actually
%   be zero, even though it should be a negative value (since A is stable).
%
%   eigsPlus solves the aboves problems and provides some additional
%   features for increased functionality and convenience.
% 
%   USAGE: 
%       d           = eigsPlus(...)
%       [V,d]       = eigsPlus(...)
%       [V,d,iters] = eigsPlus(...)
%
%   INPUT: same input arguments as eigs. 
%
%   OUTPUT: modified versions of eigs output arguments.
%
%       d       Column vector of p converged eigenvalues.  Can be empty if
%               no eigenvalues converged.  eigsPlus always returns d as a 
%               column vector, regardless of whether or not the
%               eigenvectors are also requested.
% 
%       V       Matrix of p eigenvectors corresponding to the p converged
%               eigenvalues given in d.  Can be empty if no eigenvalues
%               converged.
% 
%       iters   The number of ARPACK iterations incurred.
%
%   See also eigs, makeEigsPlus, cleanEigsPlus.
%
%
%   For comments/bug reports, please visit the eigsPlus GitLab webpage:
%   https://gitlab.com/timmitchell/eigsPlus
%
%   eigsPlus Version 2.1, 2016-2018, see AGPL license info below.
%   eigsPlus.m introduced in eigsPlus Version 2.0.
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

    persistent eigs_fn;
    
    if isempty(eigs_fn)
        if isOctave()
            eigs_fn = @wrapperOctave;
        else
            eigs_fn = getEigsPlusFunction();
        end
    end
      
    [varargout{1:nargout}] = eigs_fn(varargin{:});
end