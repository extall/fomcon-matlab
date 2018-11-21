function replaceParents(pop)

    D = pop.pop_data;
    O = pop.options;

    
    % rename for clarity
    new_fits = D.function_values_offspring;
    new_inds = D.offspring_population;

    % operation depends on the algorithm again
    switch upper(pop.algorithm)

        % Multistart algorithm
        case 'MS'
            % MS just replaces everything
            D.parent_population      = new_inds;
            D.function_values_parent = new_fits;

        % Differential Evolution
        case 'DE' % Differential Evolution
            % DE and GA both use simple greedy replacement
            better_inds                              = new_fits < pop.fitnesses;
            D.parent_population(better_inds, :)      = new_inds(better_inds, :);
            D.function_values_parent(better_inds, :) = new_fits(better_inds, :);
            
        % Genetic Algorithm
        case 'GA'  % Genetic Algorithm
            % DE and GA both use simple greedy replacement
            better_inds                              = new_fits < pop.fitnesses;
            D.parent_population(better_inds, :)      = new_inds(better_inds, :);
            D.function_values_parent(better_inds, :) = new_fits(better_inds, :);

        % Particle Swarm Optimization
        case 'PSO' % Particle Swarm Optimization

            % PSO simply replaces all parents
            D.parent_population      = new_inds;
            D.function_values_parent = new_fits;

            % update the neighbor bests
            % (this implementation is fast, not intuitive)
            % add one NaN to the new_fits array
            new_fits  = [new_fits; NaN];
            % copy neighbors
            neighbors = D.neighbors;
            % let those that are zero refer to the NaN entry
            neighbors(neighbors == 0) = size(new_fits,1); 
            % find the best ones
            [neighbor_best, ind] = min(new_fits(neighbors),[],2);
            % find those that are better
            better_neighbors = neighbor_best <  D.neighbor_best_fits;
            % no better ones might be found
            if any(better_neighbors)
                % insert function values
                D.neighbor_best_fits(better_neighbors) = ...
                    neighbor_best(better_neighbors);
                % insert individuals
                for i = (find(better_neighbors)).'
                    D.neighbor_best_inds(i, :) = ...
                        new_inds(neighbors(i, ind(i)), :);
                end
            end
            % chop off additional NaN-entry again
            new_fits = new_fits(1:end-1);

            % update the local bests
            new_locals = new_fits < D.local_best_fits;
            D.local_best_fits(new_locals, 1) = new_fits(new_locals, 1);
            D.local_best_inds(new_locals, :) = new_inds(new_locals, :);

            % update the global best
            if (min(new_fits) < D.global_best_fit)
                [D.global_best_fit, ind] = min(new_fits);
                D.global_best_ind = new_inds(ind, :);
            end

            % create random matrices
            r1  = rand(pop.size, 1);  r1 = r1(:, ones(1,pop.dimensions));
            r2  = rand(pop.size, 1);  r2 = r2(:, ones(1,pop.dimensions));
            r3  = rand(pop.size, 1);  r3 = r3(:, ones(1,pop.dimensions));

            % update velocities
            D.velocities = ...
                O.PSO.omega  *D.velocities + ...
                O.PSO.eta1   *r1.*(D.neighbor_best_inds - new_inds)+...
                O.PSO.eta2   *r2.*...
                                  bsxfun(@minus, D.global_best_ind, new_inds)+...
                O.PSO.eta3   *r3.*(D.local_best_inds - new_inds);
            
            % check the bounds
            pop.pop_data = D;
            pop.honorBounds([]);

        % Adaptive Simulated Annealing
        case 'ASA' % Adaptive Simulated Annealing

            % rename some stuff
            T       = D.temperature;
            T0      = O.ASA.T0;
            nrg     = D.function_values_offspring;
            prevnrg = D.function_values_parent;
            cool    = O.ASA.CoolingSchedule;
            iters   = pop.iterations - D.iters;

            % reject or accept the new population, according to
            % the probabilistic rule
            nrgdiff  = (prevnrg - nrg);                % energy difference
            nrgdiff  = nrgdiff / max(abs(nrgdiff(:))); % rescale the differences
            ind      = nrgdiff > 0;                    % always accept better ones
            probind  = ~ind & rand(pop.size, 1) < exp( nrgdiff/T );
                                                       % accept worse ones based on
                                                       % probabalistic rule
            swapinds = ind | probind;                  % indices to be swapped

            % apply cooling schedule
            D.temperature = max(eps,cool(T, T0, iters));

            % replace the individuals
            D.parent_population(swapinds, :) = new_inds(swapinds, :);

            % also replace the function values
            D.function_values_parent(swapinds, :) = new_fits(swapinds, :);

    end

    % copy individuals and fitnesses to respective properties
    pop.fitnesses   = D.function_values_parent;
    pop.individuals = D.parent_population;
    pop.pop_data    = D;

end











