% Helper function for recursively returning all fields in a structure
% USAGE: PATHS = GETFIELDS(STRCT)
% where STRCT is the structure you want to parse.
% NB! Returns a cell array of cell arrays, each containing a unique
% struct path. DOES NOT return the actual values at each endpoint.
%
% Added for FOMCON 1.5 update.
% A.T. 2021-12-23
function fpaths1 = getfields(strct, field_path, fpaths)
    
    % The actual path
    if nargin < 2 || isempty(field_path)
        field_path = {};
    end
    
    % All fpaths
    if nargin < 3 || isempty(fpaths)
        fpaths = {};
    end
    
    if isempty(field_path)
        fields = fieldnames(strct); % Takes into account base names
    else
        new_field = getfield(strct, field_path{:});
        if ~isstruct(new_field)
            fpaths1 = fpaths;
            fpaths1{end+1} = field_path;
            return;
        else
            fields = fieldnames(new_field);
        end
    end
        
    fpaths1 = fpaths;
    for q=1:length(fields)
        field_path1 = field_path;
        field_path1{end+1} = fields{q};
        fpaths1 = getfields(strct, field_path1, fpaths1);
    end
    
end