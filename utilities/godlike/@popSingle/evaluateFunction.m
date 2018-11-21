function evaluateFunction(pop)

    % multi-objective optimization overloads this function for
    % some initialization, but returns here later.

    % find evaluation sites
    if isempty(pop.pop_data.function_values_offspring)
        sites = 1:pop.size; % only in pop-initialization
        fvs = zeros(length(sites), 1);  % Single-objective only
    else
        sites = ~isfinite(pop.pop_data.function_values_offspring(:, 1));
        fvs = zeros(nnz(sites), ...
                    size(pop.pop_data.function_values_offspring, 2));
    end

    pop_inputs = pop.pop_data.offspring_population(sites, :);
    num_pop    = nnz(sites);

    % Evaluate all functions for each population member in
    % parallel or in serial
    if pop.options.UseParallel
        parfor pop_num = 1:num_pop
            fvs(pop_num, :) = evaluate_single_function(pop, pop_inputs(pop_num, :));
        end
    else
        for pop_num = 1:num_pop
            fvs(pop_num, :) = evaluate_single_function(pop, pop_inputs(pop_num, :));
        end
    end
    pop.pop_data.function_values_offspring(sites, :) = fvs;
    
    % NOTE: count each function call as one, even though it may return
    % multiple objectives.
    pop.funevals = pop.funevals + num_pop*numel(pop.funfcn);

end

function pop_output = evaluate_single_function(pop, pop_input)

    % TODO: move this to pop_multi and keep a simpler one here. Test
    % with demo.

    % Evaluate a single iterations of function(s), so that this can be run in
    % parallel.

    if pop.options.obj_columns
        pop_output = feval(pop.funfcn{1}, pop_input);
    else
        pop_output = zeros(numel(pop.funfcn),1);
        for ii = 1:numel(pop.funfcn)
            pop_output(ii) = feval(pop.funfcn{ii}, pop_input);
        end
    end

end
