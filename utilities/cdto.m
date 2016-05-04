function cdto(fun)
%CDTO Change Current Directory to the one which contains the specified FCN.

delim      = filesep;

% Get full path for file and get the path to the directory
full_path  = which(fun);
path_parts = explode(full_path, delim);
fold_path  = implode(path_parts(1:end-1), delim);

% Change directory
cd(fold_path);

end

