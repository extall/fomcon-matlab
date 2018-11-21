function createOffspring(pop, pool, times, FE)

    % get the size of the pool
    pool_size = length(pool);

    % rename some stuff
    parent_pop = pop.individuals(pool, :);
    parent_fit = pop.fitnesses(pool, :);

    % initialize
    newpop = zeros(pop.size, pop.dimensions);           % empty new population
    newfit = NaN(pop.size, pop.options.num_objectives); % placeholder for the sites to
                                                        % evaluate the function
    % determine which algorithm to use
    type = upper(pop.algorithm);

    % generate offspring with selected algorithm
    switch type

        % Multistart
        case 'MS'  % Multistart

            % number of function evaluations still allowed
            allowed_FE = pop.options.MaxFunEvals - FE;

            % set relaxed options
            options = optimset('MaxIter', times, ...
                               'display', 'off', ...
                               'TolFun' , 10*pop.options.TolFun,...
                               'TolX'   , 10*pop.options.TolX);

            % reinitialize newpop
            newpop = parent_pop;

            % loop through the population
            for i = 1:pop.size

                % reset options
                options = optimset(options,...
                                   'MaxFunEvals', allowed_FE);

                % optimize this individual
                [newpop(i,:),...
                 newfit(i,:),...
                 ef, ...
                 output] = fminsearch(pop.funfcn{1}, ...
                                      parent_pop(i, :),...
                                      options); %#ok<ASGLU>
                                  
                % update function evaluations
                pop.funevals = pop.funevals + output.funcCount;

                % number of function evaluations still allowed
                allowed_FE = allowed_FE - pop.funevals;

                % if the maximum has been exceeded, exit
                if allowed_FE < 0

                    % first insert result into pop
                    pop.pop_data.offspring_population      = newpop;
                    pop.pop_data.function_values_offspring = newfit;

                    % replace the parents
                    pop.replace_parents;

                    % then exit
                    break;
                end
            end

        % Differential Evolution
        case 'DE'  % Differential Evolution

            % I love DE for its elegance and simplicity, and
            % yet powerful optimization qualities

            % rename some stuff
            Flb = pop.options.DE.Flb;
            Fub = pop.options.DE.Fub;
            crossconst = pop.options.DE.CrossConst;

            % Neoteric Differential Evolution
            for i = 1:pop.size
                % random indices
                base = round(rand*(pool_size-1))+1; % RANDI is slower
                d1   = round(rand*(pool_size-1))+1;
                d2   = round(rand*(pool_size-1))+1;
                % d2 may not be equal to d1
                while (d1 == d2), d2 = round(rand*(pool_size-1))+1; end
                % DE operator
                if rand < crossconst || round(rand*(pool_size-1))+1 == i
                    % DE operator when rnd < Cr
                    F = rand*(Fub-Flb) + Flb;
                    newpop(i, :) = parent_pop(base,:) +  ...
                        F*(parent_pop(d1,:) - parent_pop(d2,:));
                else
                    % insert random parent otherwise
                    rnd_ind      = round(rand*(pool_size-1))+1;
                    newpop(i, :) = parent_pop(rnd_ind, :);
                    newfit(i, :) = parent_fit(rnd_ind, :);
                end
            end % for

        % Particle Swarm Optimization
        case 'PSO' % Particle Swarm Optimization

             % first, make sure the pool is large enough
            temp_pop = zeros(pop.size, pop.dimensions);
            if numel(parent_pop) ~= pop.size*pop.dimensions
                % insert all values from the pool
                temp_pop(1:size(parent_pop,1), :) = parent_pop;
                % insert random members from [parent_pop]
                for i = size(parent_pop,1)+1:pop.size
                    randinds = round(rand*(size(parent_pop,1)-1)+1);
                    temp_pop(i, :) = parent_pop(randinds, :);
                end
                % equate the two
                parent_pop = temp_pop;
            end

            % Creating offspring with PSO is pretty simple:
            newpop = parent_pop + pop.pop_data.velocities;

            % Since the velocity can not be zero, each individual
            % changes during the creation of offspring, so the
            % function values of ALL new individuals have to be
            % re-calculated.

        % Genetic Algorithm
        case 'GA'  % Genetic Algorithm

            % Generating offspring in GA requires quite many
            % operations. Compared to DE or PSO, it's really
            % quite messy:

            % rename some stuff
            Coding     = pop.options.GA.Coding;
            MutProb    = pop.options.GA.MutationProb;
            CrossProb  = pop.options.GA.CrossProb;
            NumBits    = pop.options.GA.NumBits;
            Binary     = false;
            Real       = true;

            if strcmpi(Coding, 'Binary')
                Binary = true;
                Real   = false;
            end

            % save signs
            signs = sign(parent_pop);

            % initialize some arrays that keep track of the signs
            child_signs = zeros(2, pop.dimensions);
            new_signs   = zeros(pop.size, pop.dimensions);

            % convert to binary
            if Binary

                % find correct multiplier
                multiplier = 1;
                temp_pop   = round(parent_pop);           % initialize
                parent_pop = abs(parent_pop);             % take absolute value
                temp_pop   = abs(temp_pop);               % take absolute value
                while (max(temp_pop(:)) <= 2^(NumBits)) && ~all(temp_pop(:) == 0)
                    multiplier = multiplier*10;           % adjust multiplier
                    temp_pop   = round(parent_pop*multiplier);% convert to integers
                end
                multiplier = multiplier/10;               % correct multiplier
                temp_pop   = round(parent_pop*multiplier);% convert to integers

                % see if selected number of bits cannot represent population
                if (multiplier==0.1) && ~all(temp_pop(:)==0)
                    error([mfilename('class') ':numbits_insufficient'], ...
                          ['Maximum value in population can not be represented '
                          'by the selected number of bits. Increase ''NumBits'' ',...
                          'option, or rescale your problem.']);
                end

                % convert each column separately
                bit_representation = false(pool_size, pop.dimensions*NumBits);

                % NOTE: MATLAB's DEC2BIN() should be avoided, as it would
                % be called in a loop and is not a builtin. Moreover, the
                % output of DEC2BIN() is a string, which is pretty
                % inconvenient for the mutation operator. Therefore, do
                % the conversion manually
                for i = 1:pop.dimensions
                    % convert this column to bits
                    bits = temp_pop(:, i);
                    bits = bits*pow2(1-NumBits:0);
                    bits = floor(bits);            % oddly enough, this is much faster
                    bits = bits - fix(bits/2)*2 ;  % than using REM(bits,2)...
                    bits = logical(bits);
                    % append in output matrix
                    bit_representation(:, NumBits*(i-1)+1:NumBits*(i-1)+NumBits) = bits;
                end
                % equate the two
                parent_pop = bit_representation;
                % redefine newpop
                newpop = false(pop.size, pop.dimensions*NumBits);
                % initialize children
                children = false(2, size(parent_pop,2));
                % and define convenient conversion array
                % (starts mattering for large population sizes)
                convert_to_dec = repmat(2.^(NumBits-1:-1:0), pop.size, 1);

            % convert to array of strings in case of real-representation
            elseif Real
                % do everything in one go with INT2STR()
                % (avoid NUM2STR, as its horrifically slowin a loop)
                % (note that the signs are still included in the array)
                % TODO: PERF (Rody Oldenhuis) int2str is also not builtin...think about this more
                real_representation = int2str(abs(parent_pop)*1e18);
                % convert to array of doubles
                real_representation = real_representation - '0';
                % equate the two
                parent_pop = real_representation;
                % initialize children
                children = zeros(2, size(parent_pop, 2));
                % redefine newpop
                newpop = zeros(pop.size, size(parent_pop, 2));
            end

            % perform crossover
            for i = 1:2:pop.size-1

                % select two parents at random
                parent_inds  = round(rand(2,1)*(pool_size-1)+1);
                parents      = parent_pop(parent_inds, :);
                if Binary
                    parent_signs = signs(parent_inds, :); end

                % crossover if a random number is less than [CrossProb]
                % otherwise, just insert the two parents into the new
                % population
                if (rand < CrossProb)

                    % select random crossover point
                    crosspoint = round(rand*(size(parents,2)-1))+1;

                    % perform crossover
                    children(1, :) = [parents(1,1:crosspoint), parents(2,crosspoint+1:end)];
                    children(2, :) = [parents(2,1:crosspoint), parents(1,crosspoint+1:end)];

                    % also keep track of the signs
                    if Binary
                        index = ceil(crosspoint/NumBits);
                        child_signs(1, :) = [parent_signs(1,1:index),...
                                             parent_signs(2,index+1:end)];
                        child_signs(2, :) = [parent_signs(2,1:index),...
                                             parent_signs(1,index+1:end)];
                    end

                    % insert children
                    newpop(i:i+1, :)    = children;
                    if Binary, new_signs(i:i+1, :) = child_signs; end

                else
                    newpop(i:i+1, :)    = parents;
                    newfit(i:i+1, :)    = parent_fit(parent_inds, :);
                    if Binary, new_signs(i:i+1, :) = parent_signs; end
                end % if
            end % for

            % if the population size is an uneven number, the last entry
            % is still open. Just stick a random parent there
            if mod(pop.size,2)
                index = round(rand*(pool_size-1)+1);
                newpop(end, :)   = parent_pop(index, :);
                newfit(end, :)   = parent_fit(index, :);
                if Binary, new_signs(end,:) = signs(index,:); end
            end

            % mutation operator
            mutate = rand(pop.size, size(parent_pop,2)) < MutProb;
            % If any individual mutates, the function has to be re-evaluated
            newfit(sum(mutate,2)>0,:) = NaN;
            % Binary coding - simply flip bits
            if Binary, newpop(mutate) = ~newpop(mutate); end
            % Real coding - select a new number from [0,9]
            if Real
                % don't mutate spaces
                space_inds = newpop(mutate) == (' '-'0');
                mutate(mutate) = ~space_inds;
                % don't mutate signs
                sign_inds = newpop(mutate) == ('-'-'0');
                mutate(mutate) = ~sign_inds;
                % random new digits
                rnd_inds = round(rand(nnz(mutate),1)*9);
                % convert to strings
                newpop(mutate) = rnd_inds;
            end

            % convert back to real numbers
            if Binary

                % initialize
                temp_pop = zeros(pop.size, pop.dimensions);

                % convert columnwise again
                for ii = 1:pop.dimensions
                    % convert column to decimal representation
                    temp_pop(:, ii) = sum(convert_to_dec.*newpop(:, 1:NumBits), 2);
                    % delete entries
                    newpop(:, 1:NumBits) = [];
                end

                % divide by multiplier, and re-assign signs
                temp_pop = temp_pop/multiplier.*new_signs;
                % assign newpop
                newpop = temp_pop;

            elseif Real

                % initialize
                temp_pop = zeros(pop.size, pop.dimensions);
                % assign space character
                space = ' '-'0';
                % append one "space" to the end
                newpop(:, end+1) = space;

                % then convert back to double, column per column
                for ii = 1:pop.dimensions
                    % trim leading "spaces"
                    while all(newpop(:,1)==space), newpop = newpop(:, 2:end); end
                    % first find one that does not begin with a "space"
                    non_space = find(newpop(:,1) ~= space, 1);
                    % find indices for the next "space"
                    space_ind = find(newpop(non_space,:) == space, 1);
                    space_ind = space_ind-1;
                    % use power trick forthe conversion
                    powers = 10.^(space_ind-1:-1:0);
                    powers = powers(ones(pop.size,1),:);
                    % remove residual spaces
                    ttemp_pop = newpop(:,1:space_ind);
                    ttemp_pop(ttemp_pop == space) = 0;
                    % insert in final array
                    temp_pop(:, ii) = sum(ttemp_pop.*powers,2)/1e18;
                    % adjust newpop
                    newpop = newpop(:, space_ind+1:end);
                end

                % assign newpop
                %newpop = signs.*temp_pop;

                % CG: removed 'signs.*' because it's not needed:
                % Real representation preserves signs?
                newpop = temp_pop;
            end

        % Adpative Simulated Annealing
        case 'ASA' % Adpative Simulated Annealing

            % Generating the new population is very straightforward
            % in ASA. It only requires changing ONE of the decision
            % variables per individual, which can easily be vectorized.

            % first, make sure the pool is large enough
            temp_pop = zeros(pop.size, pop.dimensions);
            if numel(parent_pop) ~= pop.size*pop.dimensions

                % insert all values from the pool
                temp_pop(1:size(parent_pop,1), :) = parent_pop;

                % insert random members from [parent_pop]
                for ii = size(parent_pop,1)+1:pop.size
                    randinds = round(rand*(size(parent_pop,1)-1)+1);
                    temp_pop(ii, :) = parent_pop(randinds, :);
                end

                % equate the two
                parent_pop = temp_pop;
            end

            % Use Bolzmann distribution to create new points
            rands  = sqrt(pop.pop_data.temperature)*randn(pop.size, pop.dimensions);
            newpop = parent_pop + rands;
            % at least one dimension always changes, so ALL function
            % values have to be recomputed

    end % switch

    % check constraints and boundaries after
    % offspring generation
    [newpop, newfit] = pop.honorBounds(newpop, newfit);

    % insert result into pop
    pop.pop_data.offspring_population      = newpop;
    pop.pop_data.function_values_offspring = newfit;

end 

