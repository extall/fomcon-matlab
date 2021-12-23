function indata = userdirman( cmd, filename, data )
%USERDIRMAN User directory manager for storing certain data like configs
%   CMD is command. Can be 'load', 'save' for filename/data.
%   DATANAME is the desired filename; for saving, filename is split at . to
%   create a variable name to be saved. NB! DATANAME is ALWAYS associated
%   with FILENAME.
%   DATA is the data to save as .mat file.
% The function will automatically check if USERDIR/.fomcon-matlab exists,
% and if it does not, it will create this path on the user's profile dir.
%
% Added for FOMCON 1.5 update.
% A.T. 2021-12-23

    if nargin < 2
        error('Not enough input arguments');
    end
    
    % Extract data name. Data name is always associated with file name.
    dataname = explode(filename, '.');
    var = dataname{1};

    % Empty indata by default. Check if it *isempty* to
    % learn whether reading required data was successful or not
    indata = [];
    
    % Make sure the folder exists
    loc = checkstore();

    switch(lower(cmd))
        case 'load'
            
            path = [loc filesep filename];
            % Let the user know if the file doesn't exist
            if exist(path)
                indata = load(path);
                indata = indata.(var); % Extract the actual variable
            end
                
        case 'save'
            
            % Need data argument as well
            if nargin < 3
                error('Not enough input arguments');
            end
            
            % Create a name for the variable containing data to be saved
            assignhere(var, data);
            
            save([loc filesep filename], var);
            
        otherwise
            error('Wrong command supplied to USERDIRMAN');
    end
            

end

% Check that the store is actually available, if not, create it
function loc = checkstore()
    FOMCON_DIR = '.fomcon-matlab'; % This is where the main config resides
    d = [gethomedir() filesep FOMCON_DIR];
    if ~exist(d, 'dir')
        mkdir(d);
    end
    
    % Return the location as well, for convenience
    loc = d;
end

function homedir = gethomedir()
    if ispc
        homedir = getenv('USERPROFILE');
    else
        homedir = getenv('HOME');
    end
end

