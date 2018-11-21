function print_progress(loop_index,...
                        algorithm_index,...
                        generation,...
                        single,multi,...
                        algorithms,...
                        which_ones,...
                        pop,...
                        converged,...
                        options,...
                        frac_popsize,...
                        frac_iterations,...
                        output)

    % if not converged, display all relevant information
    % from the current iteration
    if ~converged
        % create counter string
        genstr = num2str(generation);
        if strcmp(genstr,'11')||strcmp(genstr,'12')||strcmp(genstr,'13')
            counter_string = 'th';
        else
            switch genstr(end)
                case '1',  counter_string = 'st';
                case '2',  counter_string = 'nd';
                case '3',  counter_string = 'rd';
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
            if single
                fprintf(1, 'Performing single-objective optimization\n');
            elseif multi
                fprintf(1, ...
                        'Performing multi-objective optimization, with %d objectives.\n',...
                        options.num_objectives);
            end

            fprintf(1, [...
                    'Total population size is %d individuals. Lower bounds on\n',...
                    'algorithm iterations is %d, upper bound is %d. Generations\n', ...
                    'lower bound is %d, upper bound is %d.\n'], ...
                    options.num_objectives,...
                    sum(options.GODLIKE.popsize),...
                    options.GODLIKE.ItersLb, options.GODLIKE.ItersUb,...
                    options.MinIters       , options.MaxIters);

        end 

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
            end 
            fprintf(1, '· · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · ·\n');
        end 

        if multi
            fprintf(1, '%3d   %6d    %10d             %10d\n', ...
                algorithm_index, pop{loop_index}.funevals, ...
                nnz(pop{loop_index}.pop_data.front_number==0),...
                nnz(pop{loop_index}.pop_data.front_number~=0));
        elseif single
            fprintf(1, '%3d   %6d    %+1.5e  %+1.5e  %+1.5e\n',...
                algorithm_index, pop{loop_index}.funevals, ...
                min(pop{loop_index}.fitnesses),std(pop{loop_index}.fitnesses),...
                output.previous_best_funcvalues(loop_index) -...
                output.best_funcvalues(loop_index));
        end

    % if we do have convergence, just display the output message
    else
        fprintf(1, '\n');
        fprintf(1, output.message);
        fprintf(1, '\n\n');
    end 

end
















