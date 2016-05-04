function state = show_pbar(cmd)
%PBAR_SHOW Enables or disables progress bar for several computation-intense functions
%
% Usage:
%          PBAR_SHOW ON - turns progress bar on
%          PBAR_SHOW OFF - turns progress bar off
	if nargin ~= 1
		cmd = '';
	end
	
	switch cmd
		case 'on'
			if ~varexists('pbar_show_')
				assignin('base', 'pbar_show_', true);
			end
		case 'off'
			if varexists('pbar_show_')
				evalin('base', 'clear pbar_show_');
			end
        otherwise
            if varexists('pbar_show_') && evalin('base', 'pbar_show_')
                if nargout > 0
                    state = true;
                else
                    disp('Progress bar is enabled');
                end
            else
                if nargout > 0
                    state = false;
                else
                    disp('Progress bar is disabled');    
                end
            end
	end
		
end