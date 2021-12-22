function [fexist, f] = cfieldexists(structure, structure_fields)
%CFIELDEXISTS Check if a field (given by a string) in a structure exists
% ver 2.0: possible to use cell array as second argument. Should look like
%          {'My', 'Deep', 'Structure'}

    fexist = 1;  % Unless proven otherwise
    if ~iscell(structure_fields)
        fields_to_check = explode(structure_fields, '.');
    else
        fields_to_check = structure_fields;
    end
  
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

