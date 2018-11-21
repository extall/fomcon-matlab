classdef popSingle < handle
% =insert documentation here=


% Please report bugs and inquiries to:
%
% Name    : Rody P.S. Oldenhuis
% E-mail  : oldenhuis@gmail.com
% Licence : 2-clause BSD (See License.txt)

% If you find this work useful, please consider a donation:
% https://www.paypal.me/RodyO/3.5

    %% Properties

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
        eq_indices         % those indices of LU and UB, for which their values are equal
        eq_values          % the corresponding equal values
        trans_individuals  % transformed individuals
                           % (sine transformed and without the equal
                           % values)
        % contents of pop_data for single-objective optimization:
        % .parent_population
        % .offspring_population
        % .function_values_parent
        % .function_values_offspring
        % .unpenalized_function_values_parent     (for constrained optimization)
        % .unpenalized_function_values_offspring  (for constrained optimization)
        % .constraint_violations_parent           (for constrained optimization)
        % .constraint_violations_offspring        (for constrained optimization)
    end

    %% Methods
    
    % Class basics
    methods (Access = public)

        % Constructor
        function pop = popSingle(varargin)
            % (delegate to private method; improves file organisation)
            pop = pop.constructPop(varargin{:});
        end
        
    end

    methods (Access = private, Hidden)
        % constructor offloaded in separate file
        pop = constructPop(pop, varargin);
    end
    
    % Public functionality
    methods        
        % a single iteration
        iterate(pop, times, FE);
    end
    
    % Internal methods
    methods (Access = protected, Hidden)

        % wrapper function which includes the equal-valued, and
        % applies the sine-transformation
        transformed_individuals = wrapperFcn(pop, input_population, sites);

        % compute penalties and insert constraint violations in pop_data
        % (used for constrained optimizations)
        funvals = penalize(pop, funvals, convals, sites);

        % tournament selection (only for GA)
        pool = tournamentSelection(pop, pool_size, tournament_size);

        % generate new generation
        createOffspring(pop, pool, times, FE);

        % selectively replace the parent population with members
        % from the offspring population (single-objective optimization)
        replaceParents(pop);

        % evaluate the objective function(s) correctly
        evaluateFunction(pop);

        % check boundaries
        [newpop, newfit] = honorBounds(pop, newpop, newfit);

        % initialize algorithms
        initializeAlgorithms(pop);

    end 

end








