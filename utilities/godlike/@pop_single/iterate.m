function iterate(pop, times, FE)
% NOTE: [times] and [FE] are only used for the MultiStart algorithm

    % Select proper candiadates
    if strcmpi(pop.algorithm, 'GA')
        % binary tournament selection for GA
        pool = pop.tournament_selection(pop.size, 2);
    else
        % whole population otherwise
        pool = 1:pop.size;
    end

    % create offspring
    if nargin == 1
        pop.create_offspring(pool);
    else
        pop.create_offspring(pool, times, FE);
    end

    % if the algorithm is MS, this is the only step
    if strcmpi(pop.algorithm, 'MS')
        % adjust iterations
        pop.iterations = pop.iterations + times;
        return
    end

    % carefully evaluate objective function(s)
    try
        pop.evaluate_function;
    catch userFcn_ME
        pop_ME = MException('pop_single:function_doesnt_evaluate',...
                            'GODLIKE cannot continue: failure during function evaluation.');
        userFcn_ME = addCause(userFcn_ME, pop_ME);
        rethrow(userFcn_ME);
    end

    % Continue as usual
    pop.replace_parents();
    pop.iterations = pop.iterations + 1;

end % function (single iteration)