% wrapper function which includes the equal-valued, and 
% applies the sine-transformation
function transformed_individuals = wrapperFcn(pop,...
                                              input_population,...
                                              sites)
    
    % some initializations
    num_sites = nnz(sites);
    transformed_individuals = zeros(num_sites, prod(pop.orig_size));    
    eqind = ~pop.eq_indices(1:num_sites, :);
    
    % first unapply the sine-transformation    
    transformed_individuals(~eqind) = ...
        pop.lb(sites, :) + (pop.ub(sites,:) - pop.lb(sites,:)) .* ...
        (sin(input_population(sites, :)) + 1)/2;
    
    % then include the fixed values   
    transformed_individuals(eqind) = pop.eq_values(1:num_sites, :);
end
