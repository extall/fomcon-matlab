function replace_parents(pop)

    % rename for clarity
    new_fits = pop.pop_data.function_values_offspring;
    new_inds = pop.pop_data.offspring_population;

    % operation depends on the algorithm again
    switch upper(pop.algorithm)

        % Multistart algorithm
        case 'MS'
            % MS just replaces everything
            pop.pop_data.parent_population = new_inds;
            pop.pop_data.function_values_parent = new_fits;

        % Differential Evolution
        case 'DE' % Differential Evolution
            % DE and GA both use simple greedy replacement
            better_inds = new_fits < pop.fitnesses;
            pop.pop_data.parent_population(better_inds, :) = new_inds(better_inds, :);
            pop.pop_data.function_values_parent(better_inds, :) = new_fits(better_inds, :);
        % Genetic Algorithm
        case 'GA'  % Genetic Algorithm
            % DE and GA both use simple greedy replacement
            better_inds = new_fits < pop.fitnesses;
            pop.pop_data.parent_population(better_inds, :) = new_inds(better_inds, :);
            pop.pop_data.function_values_parent(better_inds, :) = new_fits(better_inds, :);

        % Particle Swarm Optimization
        case 'PSO' % Particle Swarm Optimization

            % PSO simply replaces all parents
            pop.pop_data.parent_population = new_inds;
            pop.pop_data.function_values_parent = new_fits;

            % update the neighbor bests
            % (this implementation is fast, not intuitive)
            % add one NaN to the new_fits array
            new_fits  = [new_fits; NaN];
            % copy neighbors
            neighbors = pop.pop_data.neighbors;
            % let those that are zero refer to the NaN entry
            neighbors(neighbors == 0) = size(new_fits,1); 
            % find the best ones
            [neighbor_best, ind] = min(new_fits(neighbors),[],2);
            % find those that are better
            better_neighbors = neighbor_best <  pop.pop_data.neighbor_best_fits;
            % no better ones might be found
            if any(better_neighbors)
                % insert function values
                pop.pop_data.neighbor_best_fits(better_neighbors) = ...
                    neighbor_best(better_neighbors);
                % insert individuals
                for i = (find(better_neighbors)).'
                    pop.pop_data.neighbor_best_inds(i, :) = ...
                        new_inds(neighbors(i, ind(i)), :);
                end
            end
            % chop off additional NaN-entry again
            new_fits = new_fits(1:end-1);

            % update the local bests
            new_locals = new_fits < pop.pop_data.local_best_fits;
            pop.pop_data.local_best_fits(new_locals, 1) = new_fits(new_locals, 1);
            pop.pop_data.local_best_inds(new_locals, :) = new_inds(new_locals, :);

            % update the global best
            if (min(new_fits) < pop.pop_data.global_best_fit)
                [pop.pop_data.global_best_fit, ind] = min(new_fits);
                pop.pop_data.global_best_ind = new_inds(ind, :);
            end

            % create random matrices
            r1  = rand(pop.size, 1);  r1 = r1(:, ones(1,pop.dimensions));
            r2  = rand(pop.size, 1);  r2 = r2(:, ones(1,pop.dimensions));
            r3  = rand(pop.size, 1);  r3 = r3(:, ones(1,pop.dimensions));

            % update velocities
            pop.pop_data.velocities = ...
                pop.options.PSO.omega  *pop.pop_data.velocities + ...
                pop.options.PSO.eta1   *r1.*(pop.pop_data.neighbor_best_inds - new_inds)+...
                pop.options.PSO.eta2   *r2.*...
                                  bsxfun(@minus, pop.pop_data.global_best_ind, new_inds)+...
                pop.options.PSO.eta3   *r3.*(pop.pop_data.local_best_inds - new_inds);

           % check the bounds
           pop.honor_bounds([]);

        % Adaptive Simulated Annealing
        case 'ASA' % Adaptive Simulated Annealing

            % rename some stuff
            T       = pop.pop_data.temperature;
            T0      = pop.options.ASA.T0;
            nrg     = pop.pop_data.function_values_offspring;
            prevnrg = pop.pop_data.function_values_parent;
            cool    = pop.options.ASA.CoolingSchedule;
            iters   = pop.iterations - pop.pop_data.iters;

            % reject or accept the new population, according to
            % the probabilistic rule
            nrgdiff  = (prevnrg - nrg);                % energy difference
            nrgdiff  = nrgdiff / max(abs(nrgdiff(:))); % rescale the differences
            ind      = nrgdiff > 0;                    % always accept better ones
            probind  = ~ind & rand(pop.size, 1) < exp( nrgdiff/T );
                                                       % accept worse ones based on
                                                       % probabalistic rule
            swapinds = ind|probind;                    % indices to be swapped

            % apply cooling schedule
            pop.pop_data.temperature = max(eps,cool(T, T0, iters));

            % replace the individuals
            pop.pop_data.parent_population(swapinds, :) = new_inds(swapinds, :);

            % also replace the function values
            pop.pop_data.function_values_parent(swapinds, :) = new_fits(swapinds, :);

    end % switch

    % copy individuals and fitnesses to respective properties
    pop.fitnesses   = pop.pop_data.function_values_parent;
    pop.individuals = pop.pop_data.parent_population;

end % function
