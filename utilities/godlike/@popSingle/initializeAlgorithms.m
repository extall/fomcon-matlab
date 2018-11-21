function initializeAlgorithms(pop)

    % PSO
    switch (pop.algorithm)
        case 'PSO'
        

            % initialize velocities
            % (average velocities about 20% of [[lb]-[ub]] interval)
            pop.pop_data.velocities = randn(pop.size, pop.dimensions) .* (pop.ub-pop.lb)/5;

            % rename for clarity
            NumNeighbors = pop.options.PSO.NumNeighbors;

            % initialize neighbors
            switch lower(pop.options.PSO.NetworkTopology)

                case 'star'
                    % star topology - in each star, there is one
                    % focal particle, to which the other members
                    % of the star are connected.

                    % initialize
                    all_particles = (1:pop.size).';
                    num_stars     = floor(pop.size/NumNeighbors);
                    pop.pop_data.neighbors = zeros(pop.size, NumNeighbors-1);

                    % initialize stars
                    if num_stars ~= 0
                        [dummy, focals] = sort(rand(pop.size,1));%#ok
                        focals = all_particles(focals(1:num_stars));
                        all_particles(focals) = [];
                        % select [NumNeighbors] random & unique neighbors
                        % for each focal particle
                        for i = 1:num_stars
                            % select new neighbors
                            [dummy, inds] = sort(rand(size(all_particles,1),1)); %#ok
                            new_neighs = all_particles(inds(1:NumNeighbors-1));
                            % adjust array
                            for j = 1:NumNeighbors-1
                                all_particles(all_particles == new_neighs(j)) = [];
                            end
                            % assign new neighbors to focal particle
                            pop.pop_data.neighbors(focals(i), :) = new_neighs;
                            % assign focal particle to new neighbors
                            pop.pop_data.neighbors(new_neighs, 1) = focals(i);
                        end
                    else
                    end

                    % population might be badly scaled for selected
                    % number of stars. Correct for this
                    % TODO - it works; those particles simply have no neighbors


                case 'ring'

                    % initialize
                    all_particles = (1:pop.size).';
                    num_rings     = floor(pop.size/NumNeighbors);
                    pop.pop_data.neighbors = zeros(pop.size, 2);

                    % form the ring
                    for i = 1:num_rings
                        % randomly select [NumNeighbors] particles
                        [dummy, inds] = sort(rand(size(all_particles,1),1)); %#ok
                        new_neighs = all_particles(inds(1:NumNeighbors));
                        % insert circularly shifted arrays
                        pop.pop_data.neighbors(new_neighs, :) = ...
                            [circshift(new_neighs,1), circshift(new_neighs,-1)];
                         % adjust array
                         for j = 1:NumNeighbors
                             all_particles(all_particles == new_neighs(j)) = [];
                         end
                    end

                    % population might be badly scaled for selected
                    % number of rings. Correct for this
                    % TODO - it works; those particles simply have no neighbors


                case 'fully_connected'
                    % fully connected swarm - all particles have
                    % ALL other particles as neighbor

                    % initialize
                    pop.pop_data.neighbors = zeros(pop.size, pop.size-1);

                    % fill the neighbors
                    for i = 1:pop.size
                        pop.pop_data.neighbors(i, :) = [1:i-1, i+1:pop.size];
                    end

            end % switch

            % find global best solution
            [global_best, index] = min(pop.fitnesses);
            pop.pop_data.global_best_ind = pop.individuals(index, :);
            pop.pop_data.global_best_fit = global_best;

            % initially, local best solutions are the function values themselves
            pop.pop_data.local_best_inds = pop.individuals;
            pop.pop_data.local_best_fits = pop.fitnesses;

            % find the neighbor best
            pop.pop_data.neighbor_best_fits = zeros(pop.size,1);
            pop.pop_data.neighbor_best_inds = zeros(pop.size, pop.dimensions);
            for i = 1:pop.size
                neighbors = pop.pop_data.neighbors(i, :);
                neighbors = neighbors(neighbors ~= 0);
                if isempty(neighbors), continue, end
                [neighbor_best, ind] = min(pop.fitnesses(neighbors));
                pop.pop_data.neighbor_best_fits(i, 1) = neighbor_best;
                pop.pop_data.neighbor_best_inds(i, :) = pop.individuals(ind, :);
            end

        % ASA
        case  'ASA'

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












