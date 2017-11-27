% overload from pop_single
function evaluate_function(pop)

    % CG: ugly hack because of pop_single constructor running
    % evaluations before pop_multi construction is complete
    if isempty(pop.num_objectives)
        num_objectives = pop.options.num_objectives;
    else
        num_objectives = pop.num_objectives;
    end

    if isempty(pop.pop_data.function_values_offspring)
        pop.pop_data.function_values_offspring = NaN(pop.size,...
                                                     num_objectives);
    end

    % call super
    evaluate_function@pop_single(pop);
end