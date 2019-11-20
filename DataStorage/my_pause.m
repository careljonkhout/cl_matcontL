function out = my_pause(command)
  persistent enabled
  if isempty(enabled)
    enabled = true;
  end
  if nargin == 1 
    switch command
      case 'on'
        out = enabled;
        enabled = true;
      case 'off'    
        out = enabled;
        enabled = false;
      case 'query'
        out = enabled;
      otherwise
        error('The command ''%s'' is not a valid command for my_pause', ...
                command);
    end
    return
  end
    
  if enabled
    input('Press enter to continue or ctrl-c to abort', 's')
  end
end

