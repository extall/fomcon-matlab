function elem = prefnumerize(current, vals)
% PREFNUMERIZE Convert given values to closest preferred series values

    % If values are given in series identifier, use that
    if isa(vals, 'char')
       prefvals = std_series_e(vals);
    else
       % Otherwise, user provides some values
       prefvals = vals;
    end
    
    curvals = scinum(current);
    seppref = scinum(prefvals);
    
    % Obtain closest matches
    newvals = [];
    for n=1:size(curvals,1)
        
        % Which exponent to use?
        if isa(vals, 'char')
           diffval  = abs(curvals(n,1)-seppref(:,1));
           [ig1, idx] = min(diffval);
           exponent = curvals(n,2); 
           newvals(end+1) = seppref(idx,1) * 10^exponent;
        else
           diffval  = abs(current(n)-vals);
           [ig1, idx] = min(diffval); 
           newvals(end+1) = vals(idx);
        end

    end
    
    % Return values
    elem = newvals;

end