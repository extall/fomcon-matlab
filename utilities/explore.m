function explore(fun)
%EXPLORE Run MS Windows EXPLORER and open directory with specified M file

% NB! This function will only work on Microsoft Windows OS at this time
delim      = filesep;
expr       = 'explorer "%FOLDER_PATH%"';

% Get full path for file and get the path to the directory
full_path  = which(fun);
path_parts = explode(full_path, delim);
fold_path  = implode(path_parts(1:end-1), delim);

% Replace in calling expression
expr = strrep(expr, '%FOLDER_PATH%', fold_path);

% Open directory in explorer
dos(expr);

end

