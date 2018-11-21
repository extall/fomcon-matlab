% non-dominated sort, and crowding distance assignment
function nonDominatedSort(pop)

    % this function sorts the population according to non-domination.
    % At the same time, it will compute the crowding distance for every
    % individual. It will store these values for all individuals in the
    % data structure [pop.pop_data].

    % NOTE: NON_DOMINATED_SORT is by far the most computationally costly
    % function of both POP_SINGLE and POP_MULTI. Any effort regarding
    % increasing efficiency should be focussed on this method.

    % determine if this is the first call or a subsequent call,
    % and create fitnesses and number of variables accordingly
    if ~isempty(pop.pop_data.function_values_offspring)
        if strcmpi(pop.algorithm, 'PSO')
            inds = [pop.pop_data.local_best_inds;
                pop.pop_data.offspring_population; ];
            fits = [pop.pop_data.local_best_fits;
                pop.pop_data.function_values_offspring];
        else
            inds = [pop.pop_data.parent_population;
                pop.pop_data.offspring_population; ];
            fits = [pop.pop_data.function_values_parent;
                pop.pop_data.function_values_offspring];
        end
        N = 2*pop.size;
    else
        inds = pop.individuals;
        fits = pop.fitnesses;
        N    = pop.size;
    end

    % pre-calculate some stuff for crowding distances
    crowding_dists = zeros(N, 1);              % initialize
    [sorted_fitnesses, indices] = sort(fits);  % sort fitnesses
    crowding_dists(indices(1, :))   = inf;     % always include boundaries
    crowding_dists(indices(end, :)) = inf;     % always include boundaries

    % count the number of solutions that dominate every other solution
    guy_is_dominated_by = inf(N, 1);
    for guy = 1:N

        % extract its fitnesses
        this_guy = fits(guy, :);

        % This is where the largest portion of the computational cost is,
        %
        %                      -> !!BY FAR!! <-
        %
        % BSXFUN() has been used to compare everything in one go, which
        % is also much faster than using a nested for. Note that strictly
        % speaking, TWO comparisons are required, since dominance is
        % defined as
        %
        % (xi <= yi) for all i, AND (xi < yi) for at least one value of i
        %
        % I coded the two comparisons for completeness, but in practice
        % two equal function values are almost never encountered, so I
        % left it commented. To use the stricter version, just uncomment
        % the second comparison (and the addition to the [dominated] array.
        less_eq   = bsxfun(@le, fits([1:guy-1, guy+1:N], :), this_guy);
        %less      = bsxfun(@lt, fits([1:guy-1, guy+1:N], :), this_guy);
        dominated = all(less_eq, 2);% & any(less, 2);

        % insert the domination count in dominatrix
        guy_is_dominated_by(guy) = nnz(dominated);

        % compute this guy's crowding distance
        for m = 1:pop.num_objectives
            % current sorting indices
            sort_inds = indices(:, m);
            % find this guy's index
            guys_ind = find(sort_inds == guy);
            % compute crowding distance
            % FIXED: (Seth Kenner) in some cases, the crowding distance turns
            % non-finite due to exactly-equal fitnesses (divide by zero). So,
            % we have to take proper precautions
            if isfinite(crowding_dists(guy))

                addition = (sorted_fitnesses(guys_ind+1, m) - sorted_fitnesses(guys_ind-1, m)) / ...
                           (sorted_fitnesses(end,m)         - sorted_fitnesses(1,m));

                if isfinite(addition)
                    crowding_dists(guy) = crowding_dists(guy) + addition; end

            else
                break
            end % if
        end % for

    end % for

    % create new population
    new_pop = zeros(pop.size, pop.dimensions);
    new_fit = zeros(pop.size, pop.num_objectives);
    fronts  = zeros(pop.size, 1);
    not_filled = pop.size;
    for i = 0:max(guy_is_dominated_by)
        % extract current front
        front = guy_is_dominated_by == i;
        % number of entries
        entries = nnz(front);
        % see if it still fits in the population
        if entries <= not_filled
            % if so, insert all individuals from this front in the population
            % and keep track of their fitnesses
            new_pop(pop.size-not_filled+1:pop.size-not_filled+entries, :) = inds(front, :);
            new_fit(pop.size-not_filled+1:pop.size-not_filled+entries, :) = fits(front, :);
            fronts (pop.size-not_filled+1:pop.size-not_filled+entries, 1) = i;
            % adjust number of entries that have not yet been filled
            not_filled = not_filled - entries;
        % if it does not fit, insert the remaining individuals based on
        % their crowding distance
        else
            break;
        end 
        
    end 

    % (for the first iteration, the WHOLE current population will fit
    % in the new population. Skip that case)
    if (N ~= pop.size)
        % sort last front encountered by crowding-distance
        front_inds = inds(front, :);       front_fits = fits(front, :);
        [sorted_front, indices] = sort(crowding_dists(front), 'descend');%#ok<ASGLU>
        % extract individuals & fitnesses in proper order
        sorted_inds = front_inds(indices, :); sorted_fits = front_fits(indices, :);
        % insert the remaining individuals in the new population
        new_pop(pop.size-not_filled+1:end, :) = sorted_inds(1:not_filled, :);
        new_fit(pop.size-not_filled+1:end, :) = sorted_fits(1:not_filled, :);
        fronts (pop.size-not_filled+1:end, 1) = i;
    end

    % insert results in data structure
    pop.pop_data.parent_population      = new_pop;
    pop.pop_data.function_values_parent = new_fit;
    pop.pop_data.front_number           = fronts;
    pop.pop_data.crowding_distance      = crowding_dists;

    % copy individuals and fitnesses to respective class properties
    pop.fitnesses   = pop.pop_data.function_values_parent;
    pop.individuals = pop.pop_data.parent_population;

    %  some algorithms need additional operations
    pop.updateAlgorithms();

end 








