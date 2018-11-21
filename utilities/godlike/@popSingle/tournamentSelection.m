function pool = tournamentSelection(pop, pool_size, tournament_size)

    % initialize mating pool
    pool = zeros(pool_size, 1);

    % total number of competitors
    rnd_inds = zeros(pool_size*tournament_size,1);

    % create random indices outside the loop (faster)
    for i = 1:floor(pool_size*tournament_size/pop.size)
        offset = pop.size*(i-1);
        [dummy, rnds] = sort(rand(pop.size,1));%#ok<ASGLU>
        rnd_inds(offset+1:min(end,offset+pop.size), :) = rnds(1:min(end,nnz(~rnd_inds)));
    end

    % fill the mating pool
    for ii = 1:pool_size

        % select [tournament_size] individuals at random
        inds = rnd_inds(1:tournament_size);
        rnd_inds = rnd_inds(tournament_size+1:end);

        % let them compete according to
        % (xj < yj) if fit(xj) < fit(yj)
        [best, ind] = min(pop.fitnesses(inds));%#ok<ASGLU>

        % insert the index of the best one in the pool
        pool(ii) = inds(ind);

    end

end
