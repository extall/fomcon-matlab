% perform one multi-objective iteration
function iterate(pop)

    % one iteration is:
    pop.non_dominated_sort;      % non-dominated sort
    if ~strcmpi(pop.algorithm, 'PSO')
        pool = ...               % binary tournament selection (half pop.size, not for PSO)
            pop.tournament_selection(round(pop.size/2), 2);
    else
        pool = 1:pop.size;       % the whole population for PSO
    end
    pop.create_offspring(pool);  % create new offspring
    pop.evaluate_function;       % evaluate objective function(s)

    % increase number of iterations made
    pop.iterations = pop.iterations + 1;

end % function (iterate)