function time_out = time_print_temp(time_in)
  persistent stored_time
  if nargin
    stored_time = time_in;
  end
  time_out = stored_time;
end

