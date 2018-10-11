function debug(varargin)
%
% debug routine, only for testing purposes
%

if nargin == 0
  error('debug needs at least one parameter');
end

%fprintf('DEBUG: ');

if nargin == 1
  fprintf(varargin{1});
else
  fprintf(varargin{1}, varargin{2:nargin});
end
