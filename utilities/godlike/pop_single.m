classdef pop_single < handle
% =insert documentation here=

%      Author : Rody P.S. Oldenhuis
% Affiliation : Delft University of Technology
%               Faculty of Aerospace Engineering
%               Dep. of Astrodynamics & Satellite Systems 
%     Contact : oldnhuis@dds.nl
%   Licensing/
%    (C) info : Frankly I don't care what you do with it, 
%               as long as I get some credit when you copy 
%               large portions of the code ^_^
    
    % all properties are public
    properties 
        algorithm          % type of optimization algorithm used 
        funfcn             % objective function(s)
        individuals        % members of the population
        fitnesses          % corresponding fitnesses        
        size               % population size
        lb                 % lower bounds
        ub                 % upper bounds
        orig_size          % original size of the input
        dimensions         % dimensions
        funevals   = 0;    % number of function evaluations made        
        iterations = 0;    % iterations so far performed
        options            % options structure (see function [set_options] for info)            
        pop_data           % structure to store intermediate data                           
        % contents for single-objective optimization:
        %      pop_data.parent_population
        %      pop_data.offspring_population
        %      pop_data.function_values_parent
        %      pop_data.function_values_offspring
    end
        
    % public methods
    methods (Access = public)
        
        % constructor
        function pop = pop_single(varargin)
            
            % default check
            error(nargchk(2, 7, nargin));
            
            % input is ( new [pop_data] structure, previous [population] object, options )
            % (subsequent call from GODLIKE)
            % = = = = = = = = = = = = = = = = = = = = = = = = = = 
            
            if (nargin == 3) 
                
                % assign new pop_data structure
                pop.pop_data = varargin{1};
                
                % simply copy previous object
                pop.funfcn     = varargin{2}.funfcn;       pop.iterations = varargin{2}.iterations;
                pop.algorithm  = varargin{2}.algorithm;    pop.lb         = varargin{2}.lb;
                pop.funevals   = varargin{2}.funevals;     pop.ub         = varargin{2}.ub;
                pop.dimensions = varargin{2}.dimensions;   pop.orig_size  = varargin{2}.orig_size;  
                
                % copy individuals and fitnesses
                pop.individuals = pop.pop_data.parent_population;
                pop.fitnesses   = pop.pop_data.function_values_parent;
                
                % size and options might have changed
                pop.size = size(pop.individuals, 1);%#ok
                pop.options = varargin{3}; 
                % replicate [ub] and [lb]
                pop.lb = repmat(pop.lb(1, :), pop.size, 1);    
                pop.ub = repmat(pop.ub(1, :), pop.size, 1);
                
                % Some algorithms need some lengthier initializing
                pop.initialize_algorithms;
                
                % return
                return
            end
            
            % input is ( funfcn, popsize, lb, ub, dimensions, options )
            % (initialization call from GODLIKE)
            % = = = = = = = = = = = = = = = = = = = = = = = = = = 
             
            % parse input
            % · · · · · · · · · · · · · · · · · · · · · · 
                        
            % assign input
            pop.funfcn  = varargin{1};   pop.ub         = varargin{4};
            pop.size    = varargin{2};   pop.orig_size  = varargin{5};
            pop.lb      = varargin{3};   pop.dimensions = varargin{6};
            pop.options = varargin{7};
            
            % cast funfcn to cell if necessary
            if ~iscell(pop.funfcn), pop.funfcn = {pop.funfcn}; end
                        
            % replicate [lb] and [ub] to facilitate things a bit
            % (and speed it up some more)
            pop.lb = repmat(pop.lb, pop.size, 1);   pop.ub = repmat(pop.ub, pop.size, 1);
            
            % set optimization algorithm
            pop.algorithm = pop.options.algorithm;
                             
            % Initialize population             
            % · · · · · · · · · · · · · · · · · · · · · · 
            
            % initialize population 
            pop.individuals = pop.lb + rand(pop.size, pop.dimensions) .* (pop.ub-pop.lb);
                                      
            % insert copy into info structure
            pop.pop_data.parent_population = pop.individuals;
                        
            % temporarily copy parents to offspring positions 
            pop.pop_data.function_values_offspring = [];
            pop.pop_data.offspring_population      = pop.individuals;
            
            % evaluate function for initial population (parents only)             
            pop.evaluate_function;   
                        
            % copy function values into fitnesses properties
            pop.fitnesses = pop.pop_data.function_values_offspring;
            pop.pop_data.function_values_parent = pop.fitnesses;
                        
            % delete entry again
            pop.pop_data.function_values_offspring = [];
            
            % some algorithms need some lengthier initializing
            pop.initialize_algorithms;
                                  
        end % function (constructor)  
        
        % single iteration
        function iterate(pop, times, FE)
            % [times] and [FE] are only used for the MultiStart algorithm
                        
            % select proper candiadates 
            if strcmpi(pop.algorithm, 'GA')
                pool = ...               % binary tournament selection for GA
                pop.tournament_selection(pop.size, 2); 
            else
                pool = 1:pop.size;       % whole population otherwise
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
                % then return                
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
            
            % replace the parents
            pop.replace_parents;         
                        
            % increase number of iterations made
            pop.iterations = pop.iterations + 1;
                        
        end % function (single iteration)              
                
    end % methods
    
    % % protected/hidden methods
    methods (Access = protected, Hidden)
        
        % tournament selection (only for GA)
        function pool = tournament_selection(pop, pool_size, tournament_size)
            
            % initialize mating pool
            pool = zeros(pool_size, 1);
            
            % total number of competitors
            rnd_inds = zeros(pool_size*tournament_size,1);
            
            % create random indices outside the loop (faster)
            for i = 1:floor(pool_size*tournament_size/pop.size)
                offset = pop.size*(i-1);
                [dummy, rnds] = sort(rand(pop.size,1));
                rnd_inds(offset+1:min(end,offset+pop.size), :) = rnds(1:min(end,nnz(~rnd_inds)));
            end
            
            % fill the mating pool
            for i = 1:pool_size    
                
                % select [tournament_size] individuals at random
                inds = rnd_inds(1:tournament_size);
                rnd_inds = rnd_inds(tournament_size+1:end);
                
                % let them compete according to
                % (xj < yj) if fit(xj) < fit(yj)
                [best, ind] = min(pop.fitnesses(inds));                
                
                % insert the index of the best one in the pool
                pool(i) = inds(ind);                
            
            end % for
            
        end % function (tournament selection)        
        
        % generate new generation
        function create_offspring(pop, pool, times, FE)
            
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
                    options = optimset(...
                        'MaxIter', times, ...
                        'display', 'off', ...
                        'TolFun' , 10*pop.options.TolFun,...
                        'TolX'   , 10*pop.options.TolX);
                    
                    % reinitialize newpop
                    newpop = parent_pop;
                                        
                    % loop through the population
                    for i = 1:pop.size                        
                        % reset options 
                        options = optimset(options, 'MaxFunEvals', allowed_FE);                        
                        % optimize this individual
                        [newpop(i, :), newfit(i, :), ef, output] =...
                            fminsearch(pop.funfcn{1}, parent_pop(i, :), options);
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
                            break                        
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
                        if rand < crossconst || round(rand*(pool_size-1))+1 == i;
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
                        temp_pop(1:size(parent_pop,1), :) = parent_pop; %#ok                        
                        % insert random members from [parent_pop]
                        for i = size(parent_pop,1)+1:pop.size %#ok                            
                            randinds = round(rand*(size(parent_pop,1)-1)+1); %#ok                            
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
                    if strcmpi(Coding, 'Binary')
                        Binary = true;  Real = false;
                    else
                        Binary = false; Real = true;
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
                        if (multiplier == 0.1) && ~all(temp_pop(:) == 0)
                            error('pop_single:numbits_insufficient', ...
                                  ['Maximum value in population can not be represented by the\n',...
                                  'selected number of bits. Increase ''NumBits'' option, or\n',...
                                  'rescale your problem.']);
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
                        children = false(2, size(parent_pop,2)); %#ok
                        % and define convenient conversion array
                        % (starts mattering for large population sizes)
                        convert_to_dec = repmat(2.^(NumBits-1:-1:0), pop.size, 1);
                                             
                    % convert to array of strings in case of real-representation
                    elseif Real
                        % do everything in one go with INT2STR()    
                        % (avoid NUM2STR, as its horrifically slowin a loop)
                        % (note that the signs are still included in the array)
                        real_representation = int2str(abs(parent_pop)*1e18);
                        % convert to array of doubles
                        real_representation = real_representation - '0';                       
                        % equate the two
                        parent_pop = real_representation;
                        % initialize children
                        children = zeros(2, size(parent_pop, 2)); %#ok
                        % redefine newpop
                        newpop = zeros(pop.size, size(parent_pop, 2)); %#ok
                    end
                                                            
                    % perform crossover
                    for i = 1:2:pop.size-1
                        
                        % select two parents at random
                        parent_inds  = round(rand(2,1)*(pool_size-1)+1);
                        parents      = parent_pop(parent_inds, :);
                        if Binary, parent_signs = signs(parent_inds, :); end
                                                
                        % crossover if a random number is less than [CrossProb]
                        % otherwise, just insert the two parents into the new 
                        % population 
                        if (rand < CrossProb)
                                                        
                            % select random crossover point
                            crosspoint = round(rand*(size(parents,2)-1))+1; %#ok
                                                                                                                
                            % perform crossover
                            children(1, :) = [parents(1,1:crosspoint),parents(2,crosspoint+1:end)];
                            children(2, :) = [parents(2,1:crosspoint),parents(1,crosspoint+1:end)];
                                                         
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
                    mutate = rand(pop.size, size(parent_pop,2)) < MutProb; %#ok
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
                        for i = 1:pop.dimensions                            
                            % convert column to decimal representation   
                            temp_pop(:, i) = sum(convert_to_dec.*newpop(:, 1:NumBits), 2);                           
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
                        for i = 1:pop.dimensions
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
                            temp_pop(:, i) = sum(ttemp_pop.*powers,2)/1e18;
                            % adjust newpop
                            newpop = newpop(:, space_ind+1:end);
                        end
                        % assign newpop
                        newpop = signs.*temp_pop;                        
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
                        temp_pop(1:size(parent_pop,1), :) = parent_pop; %#ok
                        
                        % insert random members from [parent_pop]
                        for i = size(parent_pop,1)+1:pop.size %#ok                            
                            randinds = round(rand*(size(parent_pop,1)-1)+1); %#ok                            
                            temp_pop(i, :) = parent_pop(randinds, :);
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
            [newpop, newfit] = pop.honor_bounds(newpop, newfit);
                         
            % insert result into pop
            pop.pop_data.offspring_population      = newpop;
            pop.pop_data.function_values_offspring = newfit;
            
        end % function (create offspring)
        
        % selectively replace the parent population with members 
        % from the offspring population (single-objective optimization)
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
                    neighbors(neighbors == 0) = size(new_fits,1); %#ok      
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
        
        % evaluate the objective function(s) correctly
        function evaluate_function(pop)
                   
            % NOTE: suited for both single and multi-objective optimization
            
            % find evaluation sites
            if isempty(pop.pop_data.function_values_offspring)
                sites = 1:pop.size; % only in pop-initialization
            else
                sites = ~isfinite(pop.pop_data.function_values_offspring(:, 1));
            end
            
            % first convert population to cell
            true_pop = reshape(pop.pop_data.offspring_population(sites, :).', ...
                [pop.orig_size,nnz(sites)]); 
            % NOTE: for-loop is faster than MAT2CELL
            cell_pop = cell(1,1,nnz(sites));
            for ii = 1:size(true_pop,3), cell_pop{1,1,ii} = true_pop(:,:,ii); end %#ok
            
            % then evaluate all functions with cellfun
            for ii = 1:numel(pop.funfcn)
                pop.pop_data.function_values_offspring(sites, ii) = ...
                    cellfun(pop.funfcn{ii}, cell_pop);
            end
            
            % update number of function evaluations
            pop.funevals = pop.funevals + ...
                nnz(sites)*size(pop.pop_data.function_values_offspring, 2);%#ok

        end % function 
           
        % check boundaries
        function [newpop, newfit] = honor_bounds(pop, newpop, newfit)
                                    
            % find violation sites
            outsiders1 = false; outsiders2 = false;
            if ~isempty(newpop)
                outsiders1 = newpop < pop.lb;
                outsiders2 = newpop > pop.ub;
            end
            
            % PSO requires more elaborate check
            if strcmpi(pop.algorithm, 'PSO')
                
                % rename for clarity
                velocity = pop.pop_data.velocities; % extract velocities
                velUb = (pop.ub - pop.lb)/5;        % upper bounds on velocity
                velLb = (pop.ub - pop.lb)/1e50;     % lower bounds on velocity
                          
                % bounce against bounds
                if any(outsiders1(:) | outsiders2(:))
                    newpop(outsiders1)   = newpop(outsiders1) - velocity(outsiders1);
                    newpop(outsiders2)   = newpop(outsiders2) - velocity(outsiders2);
                    velocity(outsiders1) = -velocity(outsiders1);
                    velocity(outsiders2) = -velocity(outsiders2);                    
                end
                
                % limit velocity
                Velsign    = sign(velocity);
                outsiders1 = abs(velocity) > abs(velUb);
                outsiders2 = abs(velocity) < abs(velLb);
                if any(outsiders1(:)) || any(outsiders2(:))
                    velocity(outsiders1) = Velsign(outsiders1).*velUb(outsiders1);
                    velocity(outsiders2) = Velsign(outsiders2).*velLb(outsiders2);
                end
                
                % re-insert velocity
                pop.pop_data.velocities = velocity;                
                
            % boundary violations in all other algorithms 
            % are simply reinitialized     
            else  
                reinit = pop.lb + rand(pop.size, pop.dimensions).*(pop.ub-pop.lb);
                if any(outsiders1(:) | outsiders2(:))
                    newpop(outsiders1) = reinit(outsiders1);
                    newpop(outsiders2) = reinit(outsiders2);
                    % also remove any function values
                    newfit(any(outsiders1,2), :) = NaN;
                    newfit(any(outsiders2,2), :) = NaN;
                end
            end % if            
        end % function
        
        % initialize algorithms
        function initialize_algorithms(pop)
            
            % PSO
            if strcmpi(pop.algorithm, 'PSO')
                                
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
                            [dummy, focals] = sort(rand(pop.size,1));                            
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
            elseif strcmpi(pop.algorithm, 'ASA')
                
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
                
            end% if
        end % function 
                
    end % methods (private)
end % classdef
