function timestamp = timestamp()
  
  current_time  = clock;
  
  year     = current_time(1);
  month    = current_time(2);
  day      = current_time(3);
  hours    = current_time(4);
  minutes  = current_time(5);
  seconds  = current_time(6);
  str_date = sprintf( '%d_%02d_%02d',   year,  month,   day);
  % Milliseconds are nice to have in the filename, in case a run follows another
  % within one second.
  str_time = sprintf('%0d_%02d_%02.3f', hours, minutes, seconds);
  timestamp = [str_date '_' str_time];
end