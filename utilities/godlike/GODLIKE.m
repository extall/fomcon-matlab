function varargout = GODLIKE(funfcn, popsize, lb, ub, varargin)
% GODLIKE           Global optimizer that combines the power
%                   of a Genetic Algorithm, Diffential Evolution,
%                   Particle Swarm Optimization and Adaptive
%                   Simulated Annealing algorithms.
%
% Usage:
%
% (Single-objective optimization)
%================================
%   sol = GODLIKE(obj_fun, popsize, lb, ub)
%   sol = GODLIKE(..., ub, which_ones)
%   sol = GODLIKE(..., which_ones, options)
%   sol = GODLIKE(..., which_ones, 'option', value, ...)
%
%   [sol, fval] = GODLIKE(...)
%   [sol, fval, exitflag] = GODLIKE(...)
%   [sol, fval, exitflag, output] = GODLIKE(...)
%
%
% (Multi-objective optimization)
% ==============================
%   sol = GODLIKE(obj_fun12..., popsize, lb, ub)
%   sol = GODLIKE({obj_fun1, obj_fun2,...}, popsize, lb, ub)
%   sol = GODLIKE(..., ub, which_ones, options)
%   sol = GODLIKE(..., which_ones, 'option', value, ...)
%
%   [sol, fval] = GODLIKE(...)
%   [..., fval, Pareto_front] = GODLIKE(...)
%   [..., Pareto_front, Pareto_Fvals] = GODLIKE(...)
%   [..., Pareto_Fvals, exitflag] = GODLIKE(...)
%   [..., exitflag, output] = GODLIKE(...)
%
%
% INPUT ARGUMENTS:
% ================
%
%   obj_fun     The objective function of which the global minimum
%               will be determined (function_handle). For multi-
%               objective optimization, several objective functions 
%               may be provided as a cell array of function handles, 
%               or alternatively, in a single function that returns
%               the different function values along the second 
%               dimension.
%               Objective functions must accept either a [popsize x
%               dimensions] matrix argument, or a [1 x dimensions] 
%               vector argument, and return a [popsize x number of 
%               objectives] matrix or [1 x number of objective] 
%               vector of associated function values (number of 
%               objectives may be 1). With the first format, the 
%               function is evaluated vectorized, in  the second 
%               case CELLFUN() is used, which is a bit slower in 
%               general.
%
%   popsize     positive integer. Indicates the TOTAL population 
%               size, that is, the number of individuals of all 
%               populations combined. 
%
%   lb, ub      The lower and upper bounds of the problem's search
%               space, for each dimension. May be scalar in case all
%               bounds in all dimensions are equal. Note that at 
%               least ONE of these must have a size of [1 x 
%               dimensions], since the problem's dimensionality is 
%               derived from it. 
%               
%   which_ones  The algorithms to be used in the optimizations. May
%               be a single string, e.g., 'DE', in which case the 
%               optimization is equal to just running a single 
%               Differential Evolution optimization. May also be a
%               cell array of strings, e.g., {'DE'; 'GA'; 'ASA'}, 
%               which uses all the indicated algorithms. When 
%               omitted or left empty, defaults to {'DE';'GA';'PSO';
%               'ASA'} (all algorithms once). 
%
%   options/    Sets the options to be used by GODLIKE. Options may
%   'option',   be a structure created by set_options, or given as 
%      value    individual ['option', value] pairs. See set_options
%               for a list of all the available options and their 
%               defaults.
%
% OUTPUT ARGUMENTS:
% =================
%
%   sol         The solution that minizes the problem globally, 
%               of size [1 x dimensions]. For multi-objective 
%               optimization, this indicates the point with the 
%               smallest distance to the (shifted) origin. 
%
%   fval        The globally minimal function value
%
%   exitflag    Additional information to facilitate fully automated
%               optimization. Negative is `bad', positive `good'. A 
%               value of '0' indicates GODLIKE did not perform any 
%               operations and exited prematurely. A value of '1' 
%               indicates normal exit conditions. A value of '-1' 
%               indicates a premature exit due to exceeding the preset
%               maximum number of function evaluations. A value of 
%               '-2' indicates that the amount of maximum GODLIKE 
%               iterations has been exceeded, and a value of '-3' 
%               indicates no optimum has been found (only for single-
%               objective optimization).
%
%   output      structure, containing much additional information 
%               about the optimization as a whole; see the manual
%               for a more detailed description. 
%
%   (For multi-objective optimization only)
%
%   Pareto_front, Pareto_Fvals
%               The full set of non-dominated solutions, and their 
%               associated function values. 
%
%   See also pop_single, pop_multi, set_options. 


%      Author : Rody P.S. Oldenhuis
% Affiliation : Delft University of Technology
%               Faculty of Aerospace Engineering
%               Dep. of Astrodynamics & Satellite Systems 
%     Contact : oldenhuis@dds.nl
%   Licensing/
%    (C) info : Frankly I don't care what you do with it, 
%               as long as I get some credit when you copy 
%               large portions of the code ^_^
    
    %% Initialize
       
    % basic check on input
    error(nargchk(4, inf, nargin));
    
    % more elaborate check on input (nested function)
    check_input;
        
    % resize and reshape boundaries and dimensions (nested function)
    [lb, ub, sze, dimensions, which_ones, options] = reformat_input(lb, ub, varargin{:});
    
    % test input objective function(s) to determine the problem's dimensions, 
    % number of objectives and proper input format (nested function)
    [options, single, multi, test_evaluations] = test_funfcn(options);           
            
    % initialize more variables
             algorithms = numel(which_ones);         % number of algorithms to use
             generation = 1;                         % this is the first generation
                    pop = cell(algorithms,1);        % cell array of [population] objects                  
    num_funevaluations  = 0;                         % number of function evaluations
    [converged, output] = check_convergence;         % initial output structure
    
    % Initially, [output] is a large structure used to move data to and from all the
    % subfunctions. Later, it is turned into the output argument [output] by removing some
    % obsolete entries from the structure. 
    
    % do an even more elaborate check (the behavior of this 
    % nested function is changed by passing the number of 
    % requested output arguments)
    check_input(nargout);
    
    %% GODLIKE loop
    
    % GODLIKE loop
    while ~converged
        
        % randomize population sizes (minimum is 5 individuals)
        frac_popsize = break_value(popsize, 5);
        
        % randomize number of iterations per algorithm
        % ([options.GODLIKE.ItersUb] is the maximum TOTAL amount 
        % of iterations that will be spent in all of the algorithms combined.)
        frac_iterations = break_value(options.GODLIKE.ItersUb, options.GODLIKE.ItersLb);
                   
        % shuffle (or initialize) populations
        pop = interchange_populations(pop);
                        
        % loop through each algorithm
        for i = 1:algorithms
            
            % perform algorithm iterations
            if strcmpi(pop{i}.algorithm, 'MS') 
            % Multi-start behaves differently; its needs to 
            % execute its iterations inside pop_single.
            
                % save previous value of number of function evaluations
                prev_FE = pop{i}.funevals;
            
                % pass data via arguments
                pop{i}.iterate(frac_iterations(i), num_funevaluations);
                
                % adjust number of function evaluations made
                num_funevaluations = num_funevaluations + pop{i}.funevals - prev_FE;
                
            else % Perform single iterations for all other algorithms 
                counter = 0; % used for single-objective optimization
                for j = 1:frac_iterations(i)
                    
                    % do single iteration on this population
                    pop{i}.iterate;
                    
                    % calculate total number of function evaluations
                    % Appareantly, pop{:}.funevals doesn't work. So
                    % we have to do a loop through all of them.
                    funevaluations = 0;
                    for k = 1:algorithms
                        if ~isempty(pop{k}),funevaluations=funevaluations+pop{k}.funevals;end
                    end % for
                    num_funevaluations = test_evaluations + funevaluations;
                    
                    % check for convergence of this iteration
                    if multi
                        % all are non-dominated, first display progress, then exit the loop
                        if all(pop{i}.pop_data.front_number == 0)
                            if ~isempty(options.display), display_progress; end, break
                        end
                    elseif single
                        % check algorithm convergence
                        [alg_converged, output, counter] = check_convergence(false,output,counter);
                        % if converged, first display progress, then exit the loop
                        if alg_converged
                            if ~isempty(options.display), display_progress; end, break
                        end
                    end % if
                    
                    % check function evaluations, and exit if it
                    % surpasses the preset maximum
                    if (num_funevaluations >= options.MaxFunEvals)
                        % also display last iteration
                        if ~isempty(options.display), display_progress; end,
                        converged = true; break,
                    end
                    
                    % display progress at every iteration
                    if ~isempty(options.display), display_progress; end
                    
                end % algorithm loop
            end
            
            % if we have convergence inside the algorithm 
            % loop, break the main loop
            if converged, break, end
            
        end % main loop
                
        % increase generation
        generation = generation + 1;
        
        % check maximum iterations
        if (generation >= options.MaxIters), converged = true; end        
        
        % check for convergence (and update output structure)
        [converged, output] = check_convergence(converged, output);
                
    end % GODLIKE loop
    
    % display final results
    if ~isempty(options.display), display_progress; end
    
    %% output values
    
    % multi-objective optimization
    if multi  
        varargout{1} = output.most_efficient_point;
        varargout{2} = output.most_efficient_fitnesses;
        varargout{3} = output.pareto_front_individuals;
        varargout{4} = output.pareto_front_fitnesses;
        varargout{5} = output.exitflag;
        % remove some fields from output structure
        output = rmfield(output, {'pareto_front_individuals','pareto_front_fitnesses',...
            'exitflag','most_efficient_point','most_efficient_fitnesses'});
        % and output what's left
        varargout{6} = output;
        
    % single-objective optimization
    elseif single 
        
        % if all went normal
        if isfield(output, 'global_best_individual')
            varargout{1} = output.global_best_individual;
            varargout{2} = output.global_best_funval;
            varargout{3} = output.exitflag;
            % remove some fields from output structure
            outpt.algorithms = output.algorithms;   outpt.funcCount = output.funcCount;
            outpt.message    = output.message;      outpt.algorithm_info = output.algorithm_info;
            outpt.iterations = output.iterations;
            % and output 
            varargout{4} = outpt;
            
        % but, no optimum might have been found
        else
            varargout{1} = NaN(1, dimensions);
            varargout{2} = NaN;
            varargout{3} = -3;
            % remove some fields from output structure
            output = rmfield(output, {'global_best_funval', 'exitflag','descent_counter',...
                'best_individuals','best_funcvalues','previous_global_best_funval',...
                'previous_best_funcvalues'});
            % adjust message
            output.message = sprintf('%s\n\n All function values encountered were INF or NaN.\n',...
                output.message);            
            % output 
            varargout{4} = output;
        end
    end
        
    %% nested functions
    
    % =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  
    % initialization shizzle
    % =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  
    
    % elaborate error trapping
    function check_input(varargin)
        if (nargin == 0)
            if isempty(funfcn)
                error('GODLIKE:function_not_defined',...
                    'GODLIKE requires at least one objective function.');
            end
            if isempty(lb)||isempty(ub)||isempty(popsize)
                error('GODLIKE:lbubpopsize_not_defined',...
                    'GODLIKE requires arguments [lb], [ub] and [popsize].');
            end
            if ~isnumeric(lb)||~isnumeric(ub)||~isnumeric(popsize)
                error('GODLIKE:lbubpopsize_not_numeric',...
                    'Arguments [lb], [ub] and [popsize] must be numeric.');
            end
            if any(~isfinite(lb)) || any(~isfinite(ub)) || ...
                    any(  ~isreal(lb)) || any(~isreal(ub))
                error('GODLIKE:lbub_not_finite',...
                    'Values for [lb] and [ub] must be real and finite.');
            end
            if ~isvector(lb) || ~isvector(ub)
                error('GODLIKE:lbub_mustbe_vector',...
                    'Arguments [lb] and [ub] must be given as vectors.');
            end
            if ~isa(funfcn, 'function_handle')                
                % might be cell array
                if iscell(funfcn)
                    for ii = 1:numel(funfcn)
                        if ~isa(funfcn{ii}, 'function_handle')
                            error('GODLIKE:funfcn_mustbe_function_handle',...
                                'All objective functions must be function handles.');
                        end
                    end                    
                % otherwise, must be function handle
                else
                    error('GODLIKE:funfcn_mustbe_function_handle',...
                        'Objective function must be given as a function handle.');
                end
            end
            if (nargin == 6) && ~isstruct(varargin{2})
                error('GODLIKE:options_mustbe_structure',...
                    'Argument [options] must be a structure.')
            end
            if any(lb > ub)
                error('GODLIKE:lb_larger_than_ub',...
                    'All entries in [lb] must be smaller than the corresponding entries in [ub].')
            end
            if ~isscalar(popsize) || ~isreal(popsize) || ~isfinite(popsize) || popsize < 0
                error('GODLIKE:popsize_is_bad',...
                    'Argument [popsize] must be a real, positive and finite scalar.')
            end
        else
            if (5*numel(which_ones) > popsize)
                error('GOLIKE:popsize_too_small',...
                    ['Each algorithm requires a population size of at least 5.\n',...
                    'Given value for [popsize] makes this impossible. Increase\n',...
                    'argument [popsize] to at least ', num2str(5*numel(which_ones)), '.']);
            end
            if (options.GODLIKE.ItersLb > options.GODLIKE.ItersUb)
                warning('GODLIKE:ItersLb_exceeds_ItersUb',...
                    ['Value of options.GODLIKE.ItersLb is larger than value of\n',...
                    'options.GODLIKE.ItersUb. Values will simply be swapped.']);
                u_b = options.GODLIKE.ItersUb;
                options.GODLIKE.ItersUb = options.GODLIKE.ItersLb;
                options.GODLIKE.ItersLb = u_b;
            end
            if (options.GODLIKE.ItersLb > options.GODLIKE.ItersUb)
                warning('GODLIKE:MaxIters_exceeds_MinIters',...
                    ['Value of options.MinIters is larger than value of\n',...
                    'options.MaxIters. Values will simply be swapped.']);
                u_b = options.MaxIters;
                options.MaxIters = options.MinIters;
                options.MinIters = u_b;
            end
            if single
                % single objective optimization has a maximum of 4 output arguments
                error(nargoutchk(0, 4, varargin{1}))
            elseif multi
                % multi-objective optimization has a maximum of 6 output arguments
                error(nargoutchk(0, 6, varargin{1}))
            end
            if strcmpi(options.display, 'plot') && single && dimensions > 2
                warning('GODLIKE:Plotting_not_possible',...
                    ['Display type was set to ''Plot'', but the number of\n',...
                    'decision variables exceeds 2. The search space can note be\n',...
                    'displayed. Set options.display to ''off'' or ''on'' to \n',...
                    '''on'' to supress this message.'])
            end
            if strcmpi(options.display, 'plot') && multi && options.num_objectives > 3
                warning('GODLIKE:Plotting_not_possible',...
                    ['Display type was set to ''Plot'', but the number of\n',...
                    'objective functions exceeds 3. The Pareto front can \n',...
                    'not be displayed. Set options.display to ''off'' or \n',...
                    '''on'' to supress this message.'])
            end
        end % if        
    end % nested function
    
    % reshape, resize and redefine input to predictable formats
    function [lb, ub, sze, dimensions, which_ones, options] = ...
             reformat_input(lb, ub, varargin)
        
        % determine which algorithms to use
        if nargin == 2 || isempty(varargin{1})% default - use all heuristic algorithms
            which_ones = {'DE';'GA';'PSO';'ASA'};
        else           % user provided algorithms - check them
            which_ones = varargin{1};
            
            % cast to cell if only one is selected
            if ischar(which_ones), which_ones = {which_ones}; end
                        
            % check the given values
            if ~iscellstr(which_ones) 
                error('GODLIKE:algorithms_must_be_character_array',...
                      'Algorithms must be given as strings in a cell-array.')
            end
            for ii = 1:numel(which_ones)
                if ~ischar(which_ones{ii})
                    error('GODLIKE:algorithm_must_be_string',...
                          ['Algorithm must be selected by one of the strings ',...
                          '''DE'', ''GA'', ''PSO'' or ''ASA''..'])
                end
                if ~strcmpi(which_ones{ii}, 'DE')  && ...
                   ~strcmpi(which_ones{ii}, 'GA')  && ...
                   ~strcmpi(which_ones{ii}, 'PSO') && ...
                   ~strcmpi(which_ones{ii}, 'ASA') && ...
                   ~strcmpi(which_ones{ii}, 'MS') 
                    error('GODLIKE:incorrect_algorithm',...
                          'Algorithm must be either ''DE'', ''GA'', ''PSO'', ''ASA'' or ''MS''.')
                end
            end
        end
        
        % set options
        if nargin <= 3, options = set_options; end   % defaults
        if nargin == 4, options = varargin{2}; end   % structure provided
        if nargin > 4 , options = set_options(varargin{2:end}); end 
                                                     % individually provided     
                                                     
        % save the original size of [lb] or [ub]        
        max_elements = max(numel(lb),numel(ub));
        if (max_elements == numel(lb)), sze = size(lb); else  sze = size(ub); end        
                                                     
        % force [lb] and [ub] to be row vectors
        lb = lb(:).';  ub = ub(:).';        
        
        % replicate one or the other when their sizes are not equal
        if ~all(size(lb) == size(ub))
            if     isscalar(lb) 
                lb = repmat(lb, size(ub));
            elseif isscalar(ub) 
                ub = repmat(ub, size(lb));
            else
                error('GODLIKE:lbub_sizes_incorrect',...
                     ['If the size of either [lb] or [ub] is equal to the problem''s dimenions\n',...
                     'the size of the other must be 1x1.'])
            end
        end
        
        % define [dimensions]
        dimensions = numel(lb);
        
    end % nested function
    
    % test the function, and determine the amount of objectives. Here 
    % it is decided whether the optimization is single-objective or 
    % multi-objective.
    function [options, single, multi, fevals] = test_funfcn(options)
                
        % initialize
        fevals = 0;
                        
        % split multi/single objective
        fun = funfcn; % make a copy
        if iscell(funfcn) && (numel(funfcn) > 1)
            % no. of objectives is simply the amount of provided objective functions
            options.num_objectives = numel(funfcn);
            % single is definitely false
            single = false;
        elseif iscell(funfcn) && (numel(funfcn) == 1)
            % single it true but might still change to false
            single = true;
            % also convert function to function_handle in this case
            funfcn = funfcn{1};
        else
            % cast fun to cell
            fun = {funfcn};
            % single is true but might still change to false
            single = true;
        end
                
        % loop through all objective functions
        % (also works for single function)
        for ii = 1:numel(funfcn)
            
            % reshape to original size
            lb_original = reshape(lb, sze);
            
            % try to evaluate the function
            try                
                % simply evaluate the function with the lower bound                
                sol = fun{ii}(lb_original);
                
                % keep track of the number of function evaluations
                fevals = fevals + 1;
                
                % see whether single must be changed to multi
                if single && (numel(sol) > 1), single = false; end
                
                % it might happen that more than one function is provided, 
                % but that one of the functions returns more than one function 
                % value. GODLIKE cannot handle that case
                if (numel(sol) > 1) && (ii > 1)
                    error('GODLIKE:multimulti_not_allowed',...
                        ['GODLIKE cannot optimize multiple multi-objective problems ',...
                        'simultaneously.\nUse GODLIKE multiple times on each of your objective ',...
                        'functions separately.\n\nThis error is generated because the first of ',...
                        'your objective functions returned\nmultiple values, while ',...
                        'you provided multiple objective functions. Only one of\nthese formats ',...
                        'can be used for multi-objective optimization, not both.'])                        
                end
                
            % if evaluating the function fails, throw error
            catch userFcn_ME                                
                pop_ME = MException('GODLIKE:function_doesnt_evaluate',...
                    'GODLIKE cannot continue: failure during function evaluation.');
                userFcn_ME = addCause(userFcn_ME, pop_ME);
                rethrow(userFcn_ME);
            end % try/catch  
            
        end % for    
        
        % see if the optimization is multi-objective
        multi = ~single;
        
    end % nested function
    
    % =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  
    % functions used in the main loop
    % =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  
    
    % break up some [value] into a vector of random integers 
    % of length [algorithms], that sums up to [value]
    function frac_value = break_value(value, Lb)
        % NOTE: The case of these variables [Lb] and [Ub] is important. 
        % The GODLIKE arguments [lb] or [ub] may get overwritten!
                
        % only one algorithm - just return value
        if algorithms == 1, frac_value = value; return; end
        
        % initially, the upper bound is the value minus
        % (algorithms-1) times the lower bound
        Ub = value - (algorithms-1)*Lb;
        
        % create array of length [algorithms] that 
        % sums to [value]        
        frac_value = zeros(algorithms, 1);
        for ii = 1:algorithms-1 % note the minus one           
            % random value (make sure it's not zero)
            rnd = 0; while (rnd == 0), rnd = round(rand*(Ub-Lb) + Lb); end
            frac_value(ii) = rnd;            
            % adjust max. value for next iteration
            Ub = round((value - sum(frac_value))/(algorithms-ii));            
        end % for
        
        % last entry is the difference of the sum of all values and the original value
        frac_value(end) = value - sum(frac_value);
        
        % sort at random
        [dummy, inds] = sort(rand(size(frac_value,1),1));
        frac_value = frac_value(inds);
        
    end % nested function
    
    % shuffle and (re)initialize the population objects
    function pop = interchange_populations(pop)
        
        % just initialize populations if this is the first iteration
        if (generation == 1)      
            for ii = 1:algorithms                
                % set algorithm for this iteration
                options.algorithm = which_ones{ii};   
                % initialize population 
                if single
                    pop{ii} = pop_single(funfcn, frac_popsize(ii), lb, ub, sze, dimensions, options);
                else
                    pop{ii} = pop_multi(funfcn, frac_popsize(ii), lb, ub, sze, dimensions, options);
                end
            end
            % we're done
            return
        end
        
        % don't shuffle if there's only one algorithm
        if (algorithms == 1), return, end
        
        % initialize
        parent_pops = zeros(popsize, dimensions);              offspring_pops = parent_pops; 
        parent_fits = zeros(popsize, options.num_objectives);  offspring_fits = parent_fits;               
        if multi
            front_numbers      = zeros(popsize, 1);  
            crowding_distances = [front_numbers;front_numbers];    
        end
        lfe1 = 0;   lfe2 = 0;    % Last Filled Entry (lfe)
        
        % extract all current populations, their function values,     
        % and other relevant information
        for ii = 1:algorithms
            
            % rename stuff for clarity
            popinfo = pop{ii}.pop_data;          popsz = pop{ii}.size;
            
            % both for single and multi-objective
            parent_pops(lfe1+1:lfe1+popsz, :)  = popinfo.parent_population;
            parent_fits(lfe1+1:lfe1+popsz, :)  = popinfo.function_values_parent;
            offspring_pops(lfe1+1:lfe1+popsz,:)= popinfo.offspring_population;
            offspring_fits(lfe1+1:lfe1+popsz,:)= popinfo.function_values_offspring;
            
            % stuff specific for multi-objective optimization
            if multi
                front_numbers(lfe1+1:lfe1+popsz, :)        = popinfo.front_number;
                crowding_distances(lfe2+1:lfe2+2*popsz, :) = popinfo.crowding_distance;
            end
            
            % update indices
            lfe1 = lfe1 + popsz;  lfe2 = lfe2 + 2*popsz;
            
        end % for
        
        % shuffle everything at random
        [dummy, rndinds] = sort(rand(popsize, 1));
        parent_pops = parent_pops(rndinds,:);    offspring_pops = offspring_pops(rndinds,:);  
        parent_fits = parent_fits(rndinds,:);    offspring_fits = offspring_fits(rndinds,:);
        if multi
            [dummy, rndinds2]  = sort(rand(2*popsize, 1));
            front_numbers      = front_numbers(rndinds,:);            
            crowding_distances = crowding_distances(rndinds2,:);
        end        
        
        % re-initialize populations accordingly
        for ii = 1:algorithms  
            
            % rename for clarity
            fp = frac_popsize(ii);
            
            % split everything up according to current [frac_popsize]
            new_popinfo.parent_population         = parent_pops(1:fp, :);
            new_popinfo.function_values_parent    = parent_fits(1:fp, :);
            new_popinfo.offspring_population      = offspring_pops(1:fp, :);
            new_popinfo.function_values_offspring = offspring_fits(1:fp, :);
            if multi
                new_popinfo.front_number      = front_numbers(1:fp, :);
                new_popinfo.crowding_distance = crowding_distances(1:fp, :);
            end % if 
            
            % change options - options for ASA are always different
            options = pop{ii}.options;   
            
            % apply re-heating
            options.ASA.T0 = options.ASA.T0 / options.ASA.ReHeating / generation;
            
            % re-initialize
            if single, pop{ii} = pop_single(new_popinfo, pop{ii}, options); 
            else       pop{ii} = pop_multi (new_popinfo, pop{ii}, options); 
            end
            
            % shrink arrays (using "... = [];" for deletion is rather slow)
            parent_pops = parent_pops(fp+1:end,:);  offspring_pops = offspring_pops(fp+1:end,:); 
            parent_fits = parent_fits(fp+1:end,:);  offspring_fits = offspring_fits(fp+1:end,:);            
            if multi
                front_numbers      = front_numbers(fp+1:end,:);
                crowding_distances = crowding_distances(fp+1:end,:);
            end % if              
        end % for        
    end % nested function
    
    % update output values, and check for convergence
    function [converged, output, counter] = ...
             check_convergence(converged, output, varargin)
        
         % some algorithms might be doubly used. 
         % save which ones they are
         persistent sames
         
         % no input - initialize
         if (nargin == 0)
             
             % initially, no convergence
             converged = false;
             
             % some algorithms might be doubly used. Find out
             % which ones, and create proper indices             
             sames = ones(algorithms, 1); 
             for ii = 1:algorithms
                 same        = strcmpi(which_ones, which_ones{ii});
                 sames(same) = 1:nnz(same); 
             end
             
             % general settings             
             output.algorithms = upper(which_ones); % algorithms used
             output.exitflag   = 0;                 % neutral exitflag
             output.message    = sprintf('No iterations have been performed.');
             output.funcCount  = 0;
             for ii = 1:algorithms
                 output.algorithm_info.(upper(which_ones{ii}))(sames(ii)).funcCount  = 0;
                 output.algorithm_info.(upper(which_ones{ii}))(sames(ii)).iterations = 0;
             end
             
             % initialize [output] for single-objective optimization
             if single
                 output.descent_counter               = 0;
                 output.global_best_individual        = NaN(1,dimensions); 
                 output.previous_global_best_individual = NaN(1,dimensions); 
                 output.global_best_funval            = inf;     
                 output.previous_global_best_funval   = inf;
                 output.best_funcvalues               = inf(1,algorithms);
                 output.previous_best_funcvalues      = inf(1,algorithms);
                 output.best_individuals              = NaN(algorithms,dimensions);
                 output.previous_best_individuals     = NaN(algorithms,dimensions);
                 for ii = 1:algorithms                      
                     output.algorithm_info.(upper(which_ones{ii}))(sames(ii)).last_population = [];
                     output.algorithm_info.(upper(which_ones{ii}))(sames(ii)).last_fitnesses  = [];                     
                 end
             end
             
             % initialize [output] for multi-objective optimization
             if multi
                 output.pareto_front_individuals = [];
                 output.pareto_front_fitnesses   = [];
                 output.most_efficient_point     = [];
                 output.most_efficient_fitnesses = [];
             end
             
             % we're done
             return
             
             % otherwise, update according to the current status of [pops]
         else
             
             % both per-algorithm and global check needs to be performed.
             % the mode of operation depends on the presence of a third
             % input argument. If given, only the current populations is
             % checked. If  omitted, all populations are checked.
             if (nargin == 3), alg_conv = true; algorithm = i; counter = varargin{1};
             else alg_conv = false;
             end
             
             % general stuff
             output.funcCount  = num_funevaluations;
             output.iterations = generation;
             for ii = 1:algorithms                 
                 output.algorithm_info.(upper(which_ones{ii}))(sames(ii)).iterations = pop{ii}.iterations;
                 output.algorithm_info.(upper(which_ones{ii}))(sames(ii)).funcCount  = pop{ii}.funevals;
             end
             
             % convergence might already have occured. Determine the reason
             if converged
                 % maximum function evaluations have been exceeded.
                 if (num_funevaluations >= options.MaxFunEvals)
                     output.exitflag = -1;
                     output.message = sprintf(['Optimization terminated:\n',...
                         ' Maximum amount of function evaluations has been reached.\n',...
                         ' Increase ''MaxFunEvals'' option.']);
                 end
                 % maximum allowable iterations have been exceeded.
                 if (generation >= options.MaxIters)
                     output.exitflag = -2;
                     output.message = sprintf(['Optimization terminated:\n',...
                         ' Maximum amount of iterations has been reached.\n',...
                         ' Increase ''MaxIters'' option.']);
                 end % if
             end % if
             
             % stuff specific for single objective optimization
             if single
                 
                 % store previous global best function value
                 output.previous_global_best_individual = output.global_best_individual;
                 output.previous_global_best_funval = output.global_best_funval;
                 output.previous_best_funcvalues    = output.best_funcvalues;
                 output.previous_best_individuals   = output.best_individuals;
                 
                 % assign global best individuals and their function
                 % values per algorithm
                 for ii = 1:algorithms
                     [output.best_funcvalues(ii), ind] = min(pop{ii}.fitnesses);
                     output.best_individuals(ii,:) = pop{ii}.individuals(ind, :);
                 end
                 % save new global best individual and function value
                 [min_func_val, index] = min(output.best_funcvalues);
                 if (output.global_best_funval > min_func_val)
                     output.global_best_funval     = min_func_val;
                     output.global_best_individual = output.best_individuals(index, :);
                 end
                 
                 % check convergence
                 if ~converged
                     
                     % per-algorithm convergence
                     if alg_conv
                         % update counter
                         if output.best_funcvalues(algorithm) < options.AchieveFunVal
                             if abs(output.previous_best_funcvalues(algorithm) - ...
                                     output.best_funcvalues(algorithm)) <= options.TolFun &&...
                                all(abs(output.previous_best_individuals(algorithm) - ...
                                     output.best_individuals(algorithm))) <= options.TolX
                                 counter = counter + 1;
                             else counter = 0;
                             end
                         end % if
                         
                         % if counter is larger than preset maximum,
                         % convergence has been achieved
                         if (counter > options.TolIters)
                             converged = true;
                         end
                         
                     % GODLIKE-convergence
                     else
                         % update counter
                         if output.global_best_funval < options.AchieveFunVal
                             if abs(output.previous_global_best_funval - ...
                                     output.global_best_funval) <= options.TolFun && ...
                                all(abs(output.previous_global_best_individual - ...
                                     output.global_best_individual)) <= options.TolX
                                 output.descent_counter = output.descent_counter + 1;
                             else output.descent_counter = 0;
                             end
                         end % if
                         
                         % if counter is larger than preset maximum, and the
                         % minimum amount of iterations has been performed,
                         % convergence has been achieved
                         if generation > options.MinIters && (output.descent_counter > 2)
                             converged = true;
                         end % if
                     end % if
                     
                     % finalize output
                     if converged && ~alg_conv
                         % insert the last population in the output
                         for ii = 1:algorithms
                             output.algorithm_info.(which_ones{ii})(sames(ii)).last_population = ...
                                 pop{ii}.individuals;
                             output.algorithm_info.(which_ones{ii})(sames(ii)).last_fitnesses = ...
                                 pop{ii}.fitnesses;
                         end
                         % finalize output structure
                         output.exitflag = 1;
                         output.message = sprintf(['Optimization terminated:\n\n',...
                             ' Coordinate differences were less than OPTIONS.TolX, and decrease\n',...
                             ' in function value was less than OPTIONS.TolFun for two consecutive\n',...
                             ' GODLIKE-iterations. GODLIKE algorithm converged without any problems.']);
                     end
                 end % if
             end % if
             
             % stuff specific for multi-objective optimization
             if multi
                 
                 % check convergence
                 if ~converged
                     % see if the minimum amount of iterations has
                     % been performed yet
                     if generation > options.MinIters
                         % test if ALL populations are non-dominated
                         all_nd = false(algorithms, 1);
                         for ii = 1:algorithms
                             all_nd(ii) = all(pop{ii}.pop_data.front_number == 0);
                         end
                         % if we have not broken prematurely, all fronts are zero, and
                         % thus we have convergence
                         if all(all_nd), converged = true; end
                         
                         % finalize output structure
                         if converged
                             % complete output structure
                             output.exitflag = 1;
                             output.message = sprintf(['Optimization terminated:\n',...
                                 'All trial solutions of all selected algorithms are non-dominated.\n',...
                                 'GODLIKE algorithm converged without any problems.']);
                         end % if (converged)
                     end % if (MinIters check)
                 end % if ~converged
                 
                 % if converged, complete output structure
                 if converged
                     % output complete Pareto front
                     for ii = 1:algorithms
                         output.pareto_front_individuals = ...
                             [output.pareto_front_individuals; pop{ii}.individuals];
                         output.pareto_front_fitnesses = ...
                             [output.pareto_front_fitnesses; pop{ii}.fitnesses];
                     end
                     % find most efficient point and its fitnesses
                     origin              = min(output.pareto_front_fitnesses);
                     shifted_fitnesses   = bsxfun(@minus, ...
                         output.pareto_front_fitnesses, origin);
                     distances_sq        = sum(shifted_fitnesses.^2,2);
                     [mindist_sq, index] = min(distances_sq);
                     output.most_efficient_point     = output.pareto_front_individuals(index, :);
                     output.most_efficient_fitnesses = output.pareto_front_fitnesses(index, :);
                 end % if converged
             end % if multi
         end % if
    end % nested function
    
    % display the algorithm's progress   
    function display_progress
        
        % if the algorithm is multistart, only print the header
                
        % current loop indices
        loop_index = i;
        %TODO - display for MS
        if ~strcmpi(pop{i}.algorithm, 'MS'), algorithm_index = j; end 
        
        % Command window
        % ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·          
        if strcmpi(options.display, 'on') || strcmpi(options.display, 'CommandWindow')
                       
            % if not converged, display all relevant information 
            % from the current iteration
            if ~converged
                % create counter string
                genstr = num2str(generation);
                if strcmp(genstr,'11')||strcmp(genstr,'12')||strcmp(genstr,'13')
                    counter_string = 'th';
                else
                    switch genstr(end)
                        case '1', counter_string = 'st';
                        case '2', counter_string = 'nd';
                        case '3', counter_string = 'rd';
                        otherwise, counter_string = 'th';
                    end
                end
                
                % display header if this is the first generation, first
                % algorithm and first algorithm iteration
                if (generation == 1) && (loop_index == 1) && (algorithm_index == 1)
                    
                    % display which algorithms
                    if algorithms == 1
                        strings = '%s population.\n';
                    elseif algorithms == 2
                        strings = '%s and %s populations.\n';
                    else   
                        strings = repmat('%s, ', 1, algorithms-1);       
                        % insert newlines if neccessary
                        if algorithms > 4
                            for ii = 16:64:length(strings)
                                strings(ii:ii+4) = '\n%s,';                                
                            end
                        end
                        strings = [strings, 'and %s populations.\n'];
                    end
                    
                    % output header
                    fprintf(1, ['\nGODLIKE optimization started with ', strings], which_ones{:});
                    
                    % display single or multi-objective optimization, and
                    % population size, iterations low and high
                    if     single
                        fprintf(1,...
                           ['Performing single-objective optimization, with total population\n'...
                            'size of %d individuals. Lower bounds on algorithm iterations\n', ...
                            'is %d, upper bound is %d.\n'], popsize, options.GODLIKE.ItersLb, ...
                            options.GODLIKE.ItersUb);
                    elseif multi
                        fprintf(1,...                            
                           ['Performing multi-objective optimization, with %d objectives.\n',...
                            'Total population size is %d individuals. Lower bounds on\nalgorithm ',...
                            'iterations is %d, upper bound is %d.\n'], options.num_objectives,...
                            popsize, options.GODLIKE.ItersLb, options.GODLIKE.ItersUb);
                    end % if
                end % if
                
                % subsequent iterations
                if (algorithm_index == 1) % display new header
                    % check if this is a new iteration
                    if (loop_index == 1)
                        fprintf(1,...
                           ['\n==============================================================\n',...
                            '                         ITERATION %d\n'],generation);
                        if single
                            fprintf(1, ...
                                '          Current global minimum: %+1.8e\n',...
                                output.global_best_funval);
                        end
                        fprintf(1, ...
                            '==============================================================\n');
                    end
                    % display new algorithm's header
                    fprintf(1,...
                       ['· · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · ·\n',...
                        '                     %s algorithm, %d%s pass\n',...
                        '               popsize: %d, max.iterations: %d\n'],...
                        which_ones{loop_index}, generation, counter_string, ...
                        frac_popsize(loop_index), frac_iterations(loop_index));
                    if multi
                        fprintf(1, ...
                            '  #  f.count      Pareto fronts       non-Pareto fronts\n');
                    elseif single
                        fprintf(1, '  #  f.count       min(F)        std(F)         descent\n');
                    end % if
                    fprintf(1, '· · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · ·\n')
                end % if
                
                if multi
                    fprintf(1, '%3d   %6d    %10d             %10d\n', ...
                        algorithm_index, pop{loop_index}.funevals, ...
                        nnz(pop{loop_index}.pop_data.front_number==0),...
                        nnz(pop{loop_index}.pop_data.front_number~=0))
                elseif single
                    fprintf(1, '%3d   %6d    %+1.5e  %+1.5e  %+1.5e\n',...
                        algorithm_index, pop{loop_index}.funevals, ...
                        min(pop{loop_index}.fitnesses),std(pop{loop_index}.fitnesses),...
                        output.previous_best_funcvalues(loop_index) -...
                        output.best_funcvalues(loop_index))
                end % if
         
            % if we do have convergence, just display the output message
            else
                fprintf(1, '\n'), fprintf(1, output.message), fprintf(1, '\n\n'),
            end % if            
        end % if (commandwindow)
                
        % Plot
        % ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  
        if strcmpi(options.display, 'Plot') 
                                
            % Check problem dimensionality (can not be larger than 2)
            if single && pop{1}.dimensions > 2, return, end
            
            % Check number of objectives (can not be larger than 3)
            if ~single && pop{1}.num_objectives > 3, return, end
            
            % initialize some stuff
            % (maximum of 16 algorithms can be displayed)
            clf, hold on, minfval = [];   
            colors = {'r.';'b.';'g.';'k.';
                      'ro';'bo';'go';'ko';
                      'rx';'bx';'gx';'kx';
                      'r^';'b^';'g^';'k^';};
                              
            % loop through all algorithms
            for ii = 1:algorithms
                
                % extract function values
                fvals = pop{ii}.fitnesses;
                                
                % single-objective
                if single
                    
                    % overall minimum and maximum function values
                    if isempty(minfval)
                        minfval = min(fvals(:));  maxfval = max(fvals(:));
                        if ~isfinite(minfval), minfval = -1e-100; end
                        if ~isfinite(maxfval), maxfval = 1e100; end
                    else
                        if min(fvals(:)) < minfval, minfval = min(fvals(:)); end
                        if max(fvals(:)) > maxfval, maxfval = max(fvals(:)); end
                    end
                    
                    % also extract individuals
                    inds = pop{ii}.individuals;
                    
                    % plot the variables versus their function value
                    if (size(inds,2) == 1)      % one dimensional
                        plot(inds, fvals, colors{ii});                        
                    elseif (size(inds,2) == 2)  % two dimensional
                        plot3(inds(:, 1), inds(:, 2), fvals, colors{ii});
                    end % if
                end % if single
                    
                % multi-objective    
                if multi                    
                    % plot the function values against each other   
                    if (size(fvals,2) == 2)     % two objectives
                        plot(fvals(:, 1), fvals(:, 2), colors{ii});
                    elseif (size(fvals,2) == 3) % three objectives
                        plot3(fvals(:, 1), fvals(:, 2), fvals(:, 3), colors{ii});
                    end % if         
                end % if
            end % for
                        
            % adjust legend entries
            legend_entries = upper(which_ones);
            if ~converged
                legend_entries{loop_index} = ...
                    [legend_entries{loop_index}, ' \bf{(evaluating)}'];
            end
            
            % make plot
            if single
                % plot & axes 
                if     (size(inds,2) == 1) % one-dimensional
                    xlabel('Decision variable x'), ylabel('Function value F(x)')
                    axis([lb(1) ub(1) minfval-1e-100 maxfval+1e-100])
                    if converged
                        plot(output.global_best_individual, output.global_best_funval, 'ko',...
                            'MarkerFaceColor', 'g', 'MarkerSize', 10)
                    end
                elseif (size(inds,2) == 2) % two-dimensional
                    xlabel('Decision variable x_1'), ylabel('Decision variable x_2')
                    zlabel('Function value F(x)'), view(30,50)
                    axis([lb(1) ub(1) lb(2) ub(2) minfval-1e-16 maxfval+1e-16])
                    if converged
                        plot3(output.global_best_individual(1),output.global_best_individual(2),...
                            output.global_best_funval, 'ko','MarkerFaceColor', 'g', ...
                            'MarkerSize', 10), view(30,50)
                    end
                end
                % make nice title
                if ~converged
                    title({'Current population versus objective function';
                        ['Generation ', num2str(generation),...
                        ', Function evaluations = ', num2str(num_funevaluations)]});
                else
                    % adjust legend
                    legend_entries{end+1} = 'global optimum';
                    % create the title
                    title({'Converged population versus objective function';
                        ['Generation ', num2str(generation),...
                        ', Function evaluations = ', num2str(num_funevaluations)];
                        '(Global optimum is the green dot)';});
                end
            elseif multi
                % plot & axes 
                if     (size(fvals,2) == 2)     % two objectives
                    xlabel('F_1(x)'), ylabel('F_2(x)')
                    if converged
                        plot(output.most_efficient_fitnesses(1),...
                            output.most_efficient_fitnesses(2), 'ko','MarkerFaceColor', 'g',...
                            'MarkerSize', 10)
                    end
                elseif (size(fvals,2) == 3)     % three objectives
                    xlabel('F_1(x)'), ylabel('F_2(x)'), zlabel('F_3(x)'), view(30,50)
                    if converged
                        plot3(output.most_efficient_fitnesses(1),...
                            output.most_efficient_fitnesses(2),output.most_efficient_fitnesses(3),...
                            'ko','MarkerFaceColor','g','MarkerSize',10), view(30,50)
                    end
                end
                % make nice title
                if ~ converged
                    title({'Current Pareto Front'; ['Generation ', num2str(generation),...
                        ', Function evaluations = ', num2str(num_funevaluations)]});
                else
                    % adjust legend
                    legend_entries{end+1} = 'most efficient';
                    % make nice title
                    title({'Final Pareto Front'; ['Generation ', num2str(generation),...
                        ', Function evaluations = ', num2str(num_funevaluations)];
                        '(Green dot is the most efficient point)'});
                end
            end
            
            % draw legend
            legend(legend_entries{:});                 
            
            % do not delay plotting
            drawnow       
                        
        end % if
    end % nested function  
    
end % function GODLIKE
