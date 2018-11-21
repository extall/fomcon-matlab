function pop = constructPop(pop, varargin)

    % nargin check
    argc = nargin - 1;        
    if verLessThan('MATLAB', '8.3')
        error(nargchk   (2,7,argc   ,'struct')); %#ok<*NCHKN>
        error(nargoutchk(0,1,nargout,'struct')); %#ok<*NCHKE>
    else
        narginchk(3,8);        
    end

    % input is ( new [pop_data] structure, previous [population] object, options )
    % (subsequent call from GODLIKE)
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    if (argc == 3)

        % assign new pop_data structure
        pop.pop_data = varargin{1};

        % simply copy previous object
        pop.funfcn      = varargin{2}.funfcn;       pop.iterations = varargin{2}.iterations;
        pop.algorithm   = varargin{2}.algorithm;    pop.lb         = varargin{2}.lb;
        pop.funevals    = varargin{2}.funevals;     pop.ub         = varargin{2}.ub;
        pop.dimensions  = varargin{2}.dimensions;   pop.orig_size  = varargin{2}.orig_size;
        %{
        pop.eq_indices  = varargin{2}.eq_indices;   pop.confcn     = varargin{2}.confcn;
        pop.constrained = varargin{2}.constrained;
        %}

        % copy individuals and fitnesses
        pop.individuals = pop.pop_data.parent_population;
        pop.fitnesses   = pop.pop_data.function_values_parent;

        % size and options might have changed
        pop.size = size(pop.individuals, 1);
        pop.options = varargin{3};
        % replicate [ub] and [lb]
        pop.lb = repmat(pop.lb(1, :), pop.size, 1);
        pop.ub = repmat(pop.ub(1, :), pop.size, 1);

        % Some algorithms need some lengthier initializing
        pop.initializeAlgorithms();
        
        return;
    end

    % input is ( funfcn, popsize, lb, ub, dimensions, options )
    % (initialization call from GODLIKE)
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    % parse input
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

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
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % initialize population
    pop.individuals = pop.lb + rand(pop.size, pop.dimensions) .* (pop.ub-pop.lb);

    % insert copy into info structure
    pop.pop_data.parent_population = pop.individuals;

    % temporarily copy parents to offspring positions
    pop.pop_data.function_values_offspring = [];
    pop.pop_data.offspring_population      = pop.individuals;

    % evaluate function for initial population (parents only)
    pop.evaluateFunction();

    % copy function values into fitnesses properties
    pop.fitnesses = pop.pop_data.function_values_offspring;
    pop.pop_data.function_values_parent = pop.fitnesses;

    % delete entry again
    pop.pop_data.function_values_offspring = [];

    % some algorithms need some lengthier initializing
    pop.initializeAlgorithms();

end
