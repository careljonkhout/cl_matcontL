% 
function printv(varargin)

if nargin==1
  fprintf('( ');
  fprintf('%f ', varargin{1});
  fprintf(')\n');
else
  fprintf('%s = ( ', varargin{1});
  fprintf('%f ', varargin{2});
  fprintf(')\n');
end


%SD:prints a vector
