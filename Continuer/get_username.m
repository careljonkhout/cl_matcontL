function username = get_username()
  [errorlevel, username] = system('echo $USER');
  if errorlevel ~= 0
    username = 'username_not_found';
    return
  end
  username = strtrim(username);
end
  
    