function updateAlgorithms(pop)

    % select which algorithm we have
    switch upper(pop.algorithm)

        % Genetic Algorithm
        case 'GA' % Genetic Algorithm
            % do nothing

        % Differential Evolution
        case 'DE' % Differential Evolution
            % do nothing

        % Partice Swarm Optimization
        case 'PSO' % Partice Swarm Optimization

            % set of current Pareto solutions
            Paretos = find(pop.pop_data.front_number == 0);

            % compute new local and global bests, compute new
            % neighbors and the new best neighbors
            for i = 1:pop.size

                % new local bests
                if all(pop.fitnesses(i, :) <= pop.pop_data.local_best_fits(i, :) ) || ...
                   any(Paretos == i)
                    % Replace local best when it's part of the new Pareto front,
                    % OR it dominates the previous local best
                    pop.pop_data.local_best_fits(i, :) = pop.fitnesses(i, :);
                    pop.pop_data.local_best_inds(i, :) = pop.individuals(i, :);
                end

                % New neighbors are [NumNeighbors] particles,
                % randomly drawn from the Pareto front
                new_neighs = Paretos(round(rand(...
                    pop.options.PSO.NumNeighbors,1)*(size(Paretos,1)-1))+1);
                pop.pop_data.neighbors(i, :) = new_neighs;

                % update the neighbor bests
                [dummy, max_dist] = max(pop.pop_data.crowding_distance(new_neighs));%#ok<ASGLU>
                best_neighbor     = new_neighs(max_dist);
                pop.pop_data.neighbor_best_fits(i, :) = ...
                    pop.fitnesses(best_neighbor, :);
                pop.pop_data.neighbor_best_inds(i, :) = ...
                    pop.individuals(best_neighbor, :);

            end % for

            % update the global best
            [maxdist, index] = max(pop.pop_data.crowding_distance(Paretos)); %#ok<ASGLU>
            if all(pop.fitnesses(Paretos(index), :) <= pop.pop_data.global_best_fit)
                % replace global best when it dominates the previous one
                pop.pop_data.global_best_fit = pop.fitnesses(Paretos(index), :);
                pop.pop_data.global_best_ind = pop.individuals(Paretos(index), :);
            end % if

        % Adaptive Simulated Annealing
        case 'ASA' % Adaptive Simulated Annealing

            % rename some stuff
            T       = pop.pop_data.temperature;
            T0      = pop.options.ASA.T0;
            cool    = pop.options.ASA.CoolingSchedule;
            iters   = pop.iterations - pop.pop_data.iters;

            % essentially, only the temperature is lowered
            pop.pop_data.temperature = max(eps, cool(T, T0, iters));

    end
    
end 














