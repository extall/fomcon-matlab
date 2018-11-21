function initializeAlgorithms(pop)

    % PSO
    switch pop.algorithm
        
        case 'PSO'

            % initialize velocities
            % (velocities are 20% of [[lb]-[ub]] interval)
            pop.pop_data.velocities = (-(pop.ub-pop.lb)/2 + ...
                rand(pop.size, pop.dimensions) .* (pop.ub-pop.lb)) / 5;

            % rename for clarity
            num_objectives = pop.options.num_objectives;
            NumNeighbors   = pop.options.PSO.NumNeighbors;

            % initialize neighbors
            % TODO - variable number of neighbors
            pop.pop_data.neighbors = round(rand(pop.size, NumNeighbors)*(pop.size-1)+1);

            % global best solution is just taken to be the first one
            pop.pop_data.global_best_ind = pop.individuals(1, :);
            pop.pop_data.global_best_fit = inf(1, num_objectives);

            % initially, local best solutions are the individuals themselves
            pop.pop_data.local_best_inds = pop.individuals;
            pop.pop_data.local_best_fits = inf(pop.size, num_objectives);

            % set the neighbor bests
            for i = 1:pop.size
                pop.pop_data.neighbor_best_fits(i, 1:num_objectives) = inf;
                pop.pop_data.neighbor_best_inds(i, :) = ...
                    pop.individuals(pop.pop_data.neighbors(i, 1), :);
            end

        case 'ASA'

            % if the initial temperature is left empty, estimate
            % an optimal one. This is simply the largest quantity
            % in (ub - lb), divided by 4, and squared. This
            % ensures that during the first few iterations,
            % particles are able to spread over the entire search
            % space; 4 = 2*2*std(randn(inf,1)).
            if isempty(pop.options.ASA.T0) && (pop.iterations == 0)
                % only do it upon initialization

                % find the maximum value
                sqrtT0dv4 = max(pop.ub(1,:) - pop.lb(1, :))/4;

                % set T0
                pop.options.ASA.T0 = sqrtT0dv4^2;

            end

            % initialize temperature
            pop.pop_data.temperature = pop.options.ASA.T0;

            % initialize iterations
            pop.pop_data.iters = pop.iterations;

    end
    
end 
