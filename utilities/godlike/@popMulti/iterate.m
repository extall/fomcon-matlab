% perform one multi-objective iteration
function iterate(pop)

    % one iteration equal:
    pop.nonDominatedSort();      % non-dominated sort
    if ~strcmpi(pop.algorithm, 'PSO')
        pool = ...               % binary tournament selection (half pop.size, not for PSO)
            pop.tournamentSelection(round(pop.size/2), 2);
    else
        pool = 1:pop.size;       % the whole population for PSO
    end
    pop.createOffspring(pool);   % create new offspring
    pop.evaluateFunction();      % evaluate objective function(s)

    
    pop.iterations = pop.iterations + 1;

end 
