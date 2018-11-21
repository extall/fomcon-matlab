
% reshape, resize and redefine input to predictable formats
function [lb,...
          ub,...
          sze,...
          popsize,...
          dimensions,...
          confcn,...
          constrained,...
          which_ones,...
          options] = reformat_input(lb,ub,...
                                    varargin)

    % set options
    if nargin <= 3, options = set_options; end                   % defaults
    if nargin == 4, options = varargin{2}; end                   % structure provided
    if nargin > 4 , options = set_options(varargin{2:end}); end  % individually provided

    % cast output functions to cell
    if isfield(options, 'OutputFcn') &&...
            isa(options.OutputFcn, 'function_handle')
        options.OutputFcn = {options.OutputFcn};
    end

    % constraint functions
    if nargin == 2 || isempty(varargin{1}) % default - no constraint function

        confcn      = {[]};
        constrained = false;

        % constraint might also be calculated inside the objective
        % function(s)
        if options.ConstraintsInObjectiveFunction > 0
            constrained = true; end

    else
        confcn      = varargin{1};
        constrained = true;

        % cast to cell if only one is selected
        if isa(confcn, 'function_handle')
            confcn = {confcn}; end

        % possible erroneous input
        if ~iscell(confcn)
            error([mfilename ':confcn_mustbe_cell_or_funhandle'], [...
                  'Constraint functions must be given as a fuction_handle, ',...
                  'or as a cell array of function_handles.']);
        end
        % FURTHER CHECKING WILL BE DONE IN TEST_FUNFCN()
    end

    % extract which algorithms to use
    which_ones = options.algorithms(:);

    % save the original size of [lb] or [ub]
    max_elements = max(numel(lb),numel(ub));
    if (max_elements == numel(lb))
        sze = size(lb);
    else
        sze = size(ub);
    end

    % force [lb] and [ub] to be row vectors
    lb = lb(:).';
    ub = ub(:).';

    % replicate one or the other when their sizes are not equal
    if ~all(size(lb) == size(ub))
        if     isscalar(lb)
            lb = repmat(lb, size(ub));
        elseif isscalar(ub)
            ub = repmat(ub, size(lb));
        else
            error([mfilename ':lbub_sizes_incorrect'], [...
                  'If the size of either [lb] or [ub] is equal to the problem''s dimenions\n',...
                  'the size of the other must be 1x1.'])
        end
    end

    % define [dimensions]
    dimensions = numel(lb);

    % total population size
    % (defaults to 25*number of dimensions)
    if isempty(options.GODLIKE.popsize)
        popsize = min(25*dimensions, 1500);
    else
        popsize = options.GODLIKE.popsize;
    end

    % check minimum popsize
    minpop = 5*numel(options.algorithms);
    if any(minpop > popsize)
        warning([mfilename ':popsize_too_small'], [...
                'Each algorithm requires a population size of at least 5.\n',...
                'Given value for [popsize] makes this impossible. Increasing\n',...
                'argument [popsize] to ', num2str(minpop), '...']);
        popsize(minpop > popsize) = minpop;
    end

    % Assign back for consistency
    options.GODLIKE.popsize = popsize;

end  
