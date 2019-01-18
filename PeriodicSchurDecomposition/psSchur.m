function [R, Q, E] = psSchur( M, Z );
%psSchur -- Periodic Schur decomposition from periodic upper Hessenberg form.
%
%   function [R, Q, E] = psSchur( M, Z )
%
%   Matlab routine input arguments:
%     * M: the array with the blocks M1 ... Mm.
%          This can be either a 2D-array with the blocks [M1 ... Mm] or a 3D array
%          with M(:,:,i) = Mi.
%          M should be in periodic upper Hessenberg form.
%     * Z: optional orthogonal transformations Z0...Z{m-1}.
%          This argument should have the same shape as M.
%          Default value: identy matrices.
%
%   Matlab output arguments:
%     * R: the array with the upper triangular factors.
%          It will have the same shape as M at exit.
%     * Q: the array with the transformation matrices.
%          It will have the same shape as M at exit.
%     * E: 3-column array.
%          + Column 1: scaled real part.
%          + Column 2: scaled imaginary part.
%          + Column 3: 10-logarithm of the scaling factor.
%          The scaling is such that the modulus of column 1 +i column 2
%          is 1, except for a zero eigenvalue. Column 3 is the 10-logarithm
%          of the modulus of the eigenvalue(s) except for zero eigenvalues.
%          If the QR iterations did not converge, a large number is returned
%          in the first and third column.
%
%   This routine computes the periodic Schur decomposition and the
%   eigenvalues of a matrix product in periodic upper Hessenberg form.
%
