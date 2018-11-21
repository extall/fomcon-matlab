% compute penalties and insert constraint violations in pop_data
% (used for constrained optimizations)
function funvals = penalize(pop, funvals, convals, sites)

    % make sure that the ones that are not violated are zero
    convals(convals <= pop.options.TolCon) = 0;

    % insert in pop.pop_data
    pop.pop_data.constraint_violations_offspring(sites, 1) = convals;
    pop.pop_data.unpenalized_function_values_offspring(sites, :) = funvals;

    % assign penalties
    if any(convals(:) > 0)
        % scaling parameter
        scale = min(1e16, 1/pop.options.TolCon);
        % replicate convals
        if (size(funvals,2) > 1)
            convals = convals(:, ones(size(funvals,2),1)); end

% TODO: (Rody Oldenhuis)
funvals = funvals + scale*convals;

%         % detect which penalties are going to overflow
%         overflows = convals > 50;
%         % penalize these with a linear penalty function
%         funvals(~overflows) = funvals(~overflows) + exp(convals(~overflows)) - 1;
%         % Use exponential penalty function for the others
%         funvals(overflows) = funvals(overflows) +  ...
%             exp(50) - 1 + scale*convals(overflows);
    end

end
