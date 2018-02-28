function plot_progress(loop_index,...
                       generation,...
                       single,multi,...
                       lb,ub,...
                       number_of_algorithms,...
                       which_ones,...
                       pop,...
                       converged,...
                       output)

    persistent ghandles
    persistent pdata
    persistent previous_legend_entries

    legend_update_needed = false;

    % Check problem dimensionality
    if single && pop{1}.dimensions     > 2  || ...
       multi  && pop{1}.num_objectives > 3
        return;
    end

    % initialize some stuff
    % (maximum of 16 algorithms can be displayed)
    if isempty(ghandles) || ...
            ~isfield(ghandles, 'fighandle') ||...
            isempty(ghandles.fighandle) || ...
            ~ishandle(ghandles.fighandle)

         % "consts"
        bgcolor = get(0, 'DefaultUiControlBackgroundColor');
        pdata = struct('number_of_algorithms', number_of_algorithms,...
                       'converged'           , converged,...
                       'lb'                  , lb,...
                       'ub'                  , ub,...
                       'pop'                 , {pop},...
                       'colors'              , {{'r.'; 'b.'; 'g.'; 'k.'
                                                 'ro'; 'bo'; 'go'; 'ko'
                                                 'rx'; 'bx'; 'gx'; 'kx'
                                                 'r^'; 'b^'; 'g^'; 'k^'}}...
        );

        ghandles.plothandles = [];
        legend_update_needed = true;

        ghandles.fighandle = figure('renderer', 'opengl',...
                                    'color'   , bgcolor,...
                                    'position', [450 150 1200 675]);
        clf, hold on
        grid on

        ghandles.gencounter_text = uicontrol('style'   , 'text',...
                                             'units'   , 'normalized',...
                                             'position', [0.02 0.96 0.20 0.03],...
                                             'string'  , ['Generation ', num2str(generation)],...
                                             'BackgroundColor'    , bgcolor,...
                                             'HorizontalAlignment', 'left');

        ghandles.fcnevals_text   = uicontrol('style'   , 'text',...
                                             'units'   , 'normalized',...
                                             'position', [0.02 0.935 0.20 0.03],...
                                             'string'  , ['Function evaluations = ', num2str(output.funcCount)],...
                                             'BackgroundColor'    , bgcolor,...
                                             'HorizontalAlignment', 'left');

        % instantiating UI objects automatically clears the toolbar; restore it
        set(ghandles.fighandle,...
            'toolbar', 'figure');

        if single
            ghandles = single_objective('init', ghandles,{}, pdata); end

        if multi
            ghandles = multi_objective('init', ghandles,{}, pdata); end

    else
        pdata.pop       = pop;
        pdata.output    = output;
        pdata.converged = converged;

        % Ready the figure for new plots
        set(0, 'currentfigure', ghandles.fighandle);
        delete(ghandles.plothandles(ishandle(ghandles.plothandles)));
        ghandles.plothandles = [];

        % Update static texts
        set(ghandles.gencounter_text, ...
            'string', ['Generation ', num2str(generation)]);

        set(ghandles.fcnevals_text,...
            'string', ['Function evaluations = ', num2str(output.funcCount)]);

    end

    % Initialize legend entries
    legend_entries = upper(which_ones);
    if ~converged
        legend_entries{loop_index} = ...
            [legend_entries{loop_index}, ' (evaluating)'];
    end

    % Plot
    if single
        [ghandles,...
         legend_entries] = single_objective('iter',...
                                            ghandles,...
                                            legend_entries,...
                                            pdata);
    elseif multi
        [ghandles,...
         legend_entries] = multi_objective('iter',...
                                           ghandles,...
                                           legend_entries,...
                                           pdata);
    end

    % Update/redraw legend
    if ~isequal(legend_entries, previous_legend_entries)
        legend_update_needed = true;  end

    if legend_update_needed
        legend([ghandles.function_outline
                ghandles.chopped_function_outline
                ghandles.plothandles(:)],...
               legend_entries{:},...
               'Location', 'NorthEast');
    end

    % Do not delay plotting
    drawnow

    % Prepare for next optimization
    previous_legend_entries = legend_entries;
    if converged
        ghandles.fighandle      = [];
        previous_legend_entries = {};
        ghandles.plothandles    = [];
    end

    single_objective('finish');
    multi_objective ('finish');

end


function [ghandles,...
         legend_entries] = single_objective(stage,...
                                            ghandles,...
                                            legend_entries,...
                                            pdata)

    persistent one_dimensional
    persistent two_dimensional
    persistent initial_draw
    persistent previous_legendentries

    switch lower(stage)
        case 'init'

            previous_legendentries = {};
            initial_draw           = true;
            legend_entries         = {};

            inds = pdata.pop{1}.individuals;

            one_dimensional = size(inds,2)==1;
            two_dimensional = size(inds,2)==2;

            % Make rough outline of function
            if one_dimensional

                x = linspace(pdata.lb(1),pdata.ub(1), 100);
                
                % NOTE: evaluate_function() cannot be used; we want to get
                % function values independent of the populations. So, we do
                % a direct evaluation of the function, so we have to be a 
                % bit careful
                f_value = pdata.pop{1}.funfcn{1}(x(:));
                assert(numel(f_value) == numel(x),...
                       [mfilename ':dimension_error'], [...
                       'Incorrect dimensions for return argument of objective ',...
                       'function.\nFor inputs of size %d-by-1, it should return ',...
                       'an array of size %d-by-1.\nHowever, it returned an array ',...
                       'of size %d-by-%d.'],...
                       numel(x),numel(x), size(f_value,1), size(f_value, 2));                
                y = reshape(f_value, size(x));

                ghandles.function_outline = plot(x,y,...
                                                 'color', [0.9 0.5 0.5]);

                ghandles.chopped_function_outline = plot(x, NaN(size(y)),...
                                                         'color', [0.9 0.9 0.9]);

            elseif two_dimensional

                [x,y] = meshgrid(linspace(pdata.lb(1),pdata.ub(1), 100),...
                                 linspace(pdata.lb(2),pdata.ub(2), 100));
                             
                % NOTE: evaluate_function() cannot be used; we want to get
                % function values independent of the populations. So, we do
                % a direct evaluation of the function, so we have to be a 
                % bit careful
                f_value = pdata.pop{1}.funfcn{1}([x(:) y(:)]);
                assert(numel(f_value) == numel(x),...
                       [mfilename ':dimension_error'], [...
                       'Incorrect dimensions for return argument of objective ',...
                       'function.\nFor inputs of size %d-by-1, it should return ',...
                       'an array of size %d-by-1.\nHowever, it returned an array ',...
                       'of size %d-by-%d.'],...
                       numel(x),numel(x), size(f_value,1), size(f_value, 2)); 
                z = reshape(f_value, size(x));
                
                ghandles.function_outline = surf(x,y,z, ...
                                                 'edgecolor', [.8 .8 .8],...
                                                 'facecolor', 'b',...
                                                 'facealpha', 0.25);

                ghandles.chopped_function_outline = surf(x,y,NaN(size(z)), ...
                                                         'edgecolor', 'none',...
                                                         'facecolor', 'r',...
                                                         'facealpha', 0.1);
            end

        case 'iter'

            % update legend entries
            legend_entries = ['Function outline'
                              'Unlikely parts'
                              legend_entries];

            update_axes = false;
            if ~isequal(previous_legendentries,legend_entries)
                previous_legendentries = legend_entries;
                update_axes = true;
            end

            % Plot all function values
            for ii = 1:pdata.number_of_algorithms

                fvals = pdata.pop{ii}.fitnesses;
                inds  = pdata.pop{ii}.individuals;

                % plot the variables versus their function value
                if one_dimensional
                    ghandles.plothandles(end+1) = plot(inds, fvals, pdata.colors{ii},...
                                                       'MarkerSize', 15);
                elseif two_dimensional
                    ghandles.plothandles(end+1) = plot3(inds(:, 1), inds(:, 2), fvals, pdata.colors{ii},...
                                                        'MarkerSize', 15);
                end % if

            end % for

            % plot & axes
            if one_dimensional

                if initial_draw
                    xlabel('x')
                    ylabel('F(x)')
                end

                if update_axes
                    M = max(cellfun(@(x)max(x.fitnesses), pdata.pop));

                    outline         = get(ghandles.function_outline, 'ydata');
                    chopped_outline = get(ghandles.chopped_function_outline, 'ydata');

                    chopped_outline(outline > M) = outline(outline > M);
                    outline(outline > M) = NaN;

                    set(ghandles.function_outline, 'ydata', outline);
                    set(ghandles.chopped_function_outline, 'ydata', chopped_outline);
                end

                if pdata.converged
                    ghandles.plothandles(end+1) = plot(pdata.output.global_best_individual,...
                                                       pdata.output.global_best_funval,...
                                                       'ko',...
                                                       'MarkerFaceColor', 'g', ...
                                                       'MarkerSize'     , 20);
                    axis tight
                end

            elseif two_dimensional

                if initial_draw
                    xlabel('x_1')
                    ylabel('x_2')
                    zlabel('F(x)')

                    view(-45,+10)
                end

                if update_axes

                    M = max(cellfun(@(x)max(x.fitnesses), pdata.pop));

                    outline         = get(ghandles.function_outline, 'zdata');
                    chopped_outline = get(ghandles.chopped_function_outline, 'zdata');

                    chopped_outline(outline > M) = outline(outline > M);
                    outline(outline > M) = NaN;

                    set(ghandles.function_outline, 'zdata', outline);
                    set(ghandles.chopped_function_outline, 'zdata', chopped_outline);
                end

                if pdata.converged
                    ghandles.plothandles(end+1) = plot3(pdata.output.global_best_individual(1),...
                                                        pdata.output.global_best_individual(2),...
                                                        pdata.output.global_best_funval,...
                                                        'ko',...
                                                        'MarkerFaceColor', 'g', ...
                                                        'MarkerSize'     , 25);
                    view(-45,+10)
                    axis tight
                end
            end

            % Update title
            if ~pdata.converged
                if initial_draw
                    title('Current population versus objective function'); end

            else
                % update legend and title
                legend_entries{end+1} = 'global optimum';
                title({'Converged population versus objective function';
                       '(Global optimum is the green dot)';});
            end

            initial_draw = false;

        case 'finish'
            return;
    end
end


function [ghandles,...
          legend_entries] = multi_objective(stage,...
                                            ghandles,...
                                            legend_entries,...
                                            pdata)

    persistent two_objectives
    persistent three_objectives
    persistent initial_draw
    persistent previous_legendentries

    switch lower(stage)

        case 'init'
            initial_draw   = true;
            legend_entries = {};

            fvals = pdata.pop{1}.fitnesses;

            two_objectives   = (size(fvals,2) == 2);
            three_objectives = (size(fvals,2) == 3);

            ghandles.function_outline         = [];
            ghandles.chopped_function_outline = [];

        case 'iter'

            % Plot objectives
            if ~pdata.converged

                for ii = 1:pdata.number_of_algorithms

                    fvals = pdata.pop{ii}.fitnesses;

                    % plot the function values against each other
                    if two_objectives
                        ghandles.plothandles(end+1) = plot(fvals(:, 1), fvals(:, 2), pdata.colors{ii},...
                                                           'MarkerSize', 15);
                    elseif three_objectives
                        ghandles.plothandles(end+1) = plot3(fvals(:, 1), fvals(:, 2), fvals(:, 3), pdata.colors{ii},...
                                                            'MarkerSize', 15);
                    end % if

                end % for
            else
                legend_entries = {'Pareto front'};

                fvals = cellfun(@(x) x.fitnesses,...
                                pdata.pop,...
                                'UniformOutput', false);
                fvals = cat(1, fvals{:});

                % plot the function values against each other
                if two_objectives
                    ghandles.plothandles(end+1) = plot(fvals(:, 1), fvals(:, 2), pdata.colors{1},...
                                                       'MarkerSize', 15);
                elseif three_objectives
                    ghandles.plothandles(end+1) = plot3(fvals(:, 1), fvals(:, 2), fvals(:, 3), pdata.colors{1},...
                                                        'MarkerSize', 15);
                end % if

            end % if

            update_axes = false;
            if ~isequal(previous_legendentries,legend_entries)
                previous_legendentries = legend_entries;
                update_axes = true;
            end

            % plot & axes
            if two_objectives

                if initial_draw
                    xlabel('F_1(x)')
                    ylabel('F_2(x)')
                end

                if update_axes
                    % TODO: (Rody Oldenhuis)
                end

                if pdata.converged
                    ghandles.plothandles(end+1) = plot(pdata.output.most_efficient_fitnesses(1),...
                                                       pdata.output.most_efficient_fitnesses(2),...
                                                       'ko',...
                                                       'MarkerFaceColor', 'g',...
                                                       'MarkerSize'     , 25);
                    axis tight
                end

            elseif three_objectives

                if initial_draw
                    xlabel('F_1(x)')
                    ylabel('F_2(x)')
                    zlabel('F_3(x)')

                    view(-45,+10)
                    axis tight
                end

                if update_axes
                    axis tight
                end

                if pdata.converged
                    ghandles.plothandles(end+1) = plot3(pdata.output.most_efficient_fitnesses(1),...
                                                        pdata.output.most_efficient_fitnesses(2),...
                                                        pdata.output.most_efficient_fitnesses(3),...
                                                        'ko',...
                                                        'MarkerFaceColor', 'g',...
                                                        'MarkerSize'     , 25);
                    view(-45,+10)
                    axis tight
                end

            end

            % Update title
            if ~pdata.converged
                if initial_draw
                    title('Current Pareto Front'); end
            else
                % update legend and title
                legend_entries{end+1} = 'most efficient';
                title({'Final Pareto Front';
                      '(Green dot is the most efficient point)'});
            end

            initial_draw = false;

        case 'finish'
            return;
    end

end
