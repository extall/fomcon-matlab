function hms_string = stohms(sec)
% STOHMS Returns a string in hh:mm:ss format from given number of seconds
            
            hour    = fix(sec/3600);
            sec     = sec - hour * 3600;
			minute  = fix(sec/60);
			sec     = sec - 60 * minute;
			second  = fix(sec);
			
			hms_string =    [num2str(hour, '%02d') ':' ...
                             num2str(minute, '%02d') ':' ...
							 num2str(second, '%02d')];

end

