function [fexist, f] = cfieldexists(structure, structure_fields)
%CFIELDEXISTS Check if a field (given by a string) in a structure exists

    fexist = 1;  % Unless proven otherwise
    fields_to_check = explode(structure_fields, '.');
    k=1;
    
    substructure = structure;
    while (fexist && k<=length(fields_to_check))
        % Not a structure or Field not found
        if ~isstruct(substructure) || ...
           ~(any(strcmp(fieldnames(substructure), fields_to_check{k})))
            fexist = 0;
        else
            substructure = substructure.(fields_to_check{k});
        end
        k=k+1;
    end
    
    % Return the field, if requested
    if fexist
        f = substructure;
    else
        f = [];
    end
    
end

