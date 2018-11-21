% break up some [value] into a vector of random integers
% of length [algorithms], that sums up to [value]
function frac_value = break_value(value,...
                                  lb,...
                                  number_of_algorithms)
        
    % only one algorithm - just return value
    if number_of_algorithms == 1
        frac_value = value;
        return;
    end

    % initially, the upper bound is the value minus
    % (algorithms-1) times the lower bound
    ub = value - (number_of_algorithms-1)*lb;

    % create array of length [algorithms] that
    % sums to [value]
    frac_value = zeros(number_of_algorithms, 1);
    for ii = 1:number_of_algorithms-1 % note the minus one

        % random value (make sure it's not zero)
        rnd = 0;
        while (rnd == 0)
            rnd = round(rand*(ub-lb) + lb); end
        frac_value(ii) = rnd;

        % adjust max. value for next iteration
        ub = round((value - sum(frac_value))/(number_of_algorithms-ii));
        
    end 

    % last entry is the difference of the sum of all values and the original value
    frac_value(end) = value - sum(frac_value);

    % sort at random
    [dummy, inds] = sort(rand(size(frac_value,1),1)); %#ok<ASGLU>
    frac_value    = frac_value(inds);

end 
