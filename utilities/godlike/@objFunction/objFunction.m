% NOTE: evalutating an objFunction() object, means you're passing it the
% TRANSFORMED X, which the class will UNTRANSFORM before evaluation. After
% convergence of the algorithm you're using it with, use
% getRealX(solution_X) to get the UNTRANSFORMED values, and use 
% obj.true_Y_at_last_evaluation to get the unpenalized function value

classdef objFunction < handle
    
    %% Properties
            
    properties (SetAccess = private)
        
        % NOTE: (Rody Oldenhuis) property validation functions can't be defined
        % on immutable properties in versions older than R2017b.

        % The actual function
        funfcn = @(X)X

        % Its constraints
        A       = []
        b       = []
        Aeq     = []
        beq     = []
        lb      = []
        ub      = []        
        intcon  = []
        nonlcon = {} 
        
        % flags & parameters
        tolfun = 1e-6
        tolcon = 1e-6
        tolx   = 1e-6
        
    end
        
    properties (SetAccess = private)
        
        % Variables
        objfcn_evaluations = 0                       
        confcn_evaluations = 0
        
        original_X_size    
        transformed_X_size 
        
        true_X_at_last_evaluation 
        true_Y_at_last_evaluation 
        
        transformed_X_at_last_evaluation 
        transformed_Y_at_last_evaluation 

        % Parameters
        nonlinear_constraints_from_objective_function = false
        number_of_objectives = 1
        
    end
    
    properties (Access = private)
        
        % Flags
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        has_been_evaluated = false
        is_constrained     = false
        
        defines_multiple_objectives = false
              
        % Indices
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        
        % Bound constraints
        nf_lb    % unconstrained towards negative infinity
        nf_ub    % unconstrained towards positive infinity
        unconst  % unconstrained in both directions
        
        fix_var  % fixed variables; lb==ub
        lb_only  % only non-finite, non-fixed constraints on LB side
        lb_ub    % non-finite, non-fixed constraints on LB and UB side
        ub_only  % only non-finite, non-fixed constraints on UB side
        
    end
    
    
    %% Methods    

    % Class basics
    methods
        
        % Constructor
        function obj = objFunction(varargin)

            
            % Initialize
            % ------------------------------------------------------------------
            
            % Allow empty objects for initialization purposes
            if nargin==0
                return; end

            % Parse the parameter/value pairs
            assert(mod(nargin,2)==0,...
                   [mfilename('class') ':pvpairs_expected'], [...
                   'Objective function object must be constructed via ',...
                   'parameter/value pairs.']);

            parameters = varargin(1:2:end);
            values     = varargin(2:2:end);

            for p = 1:numel(parameters)

                parameter = parameters{p};
                value     = values{p};

                switch lower(parameter)

                    % The objective function(s)
                    case {'function' 'fcn' 'funfcn' 'objective_function',...
                          'cost_function' 'objective' 'cost'}
                      
                        if ~iscell(value)
                            value = {value}; end
                                              
                        obj.funfcn = cellfun(@(x)verify_function_handle(x, 'objective function'),...
                                             value,...
                                             'UniformOutput', false);
                                                        
                    % Parameters
                    case 'tolfun', obj.tolfun = verify_scalar(value, 'tolfun');
                    case 'tolcon', obj.tolcon = verify_scalar(value, 'tolcon');
                    case 'tolx'  , obj.tolx   = verify_scalar(value, 'tolx'  );
                        

                    % Inequality constraints
                    % (checks are done later)
                    case 'a', obj.A = verify_matrix(value, 'A');
                    case 'b', obj.b = verify_vector(value, 'b');

                    % Equality constraints
                    % (checks are done later)
                    case {'aeq' 'a_eq' 'ae'}, obj.Aeq = verify_matrix(value, 'A_eq');
                    case {'beq' 'b_eq' 'be'}, obj.beq = verify_vector(value, 'b_eq');

                    % Bound constraints
                    % (checks are done later)
                    case {'lb' 'lower' 'lower_bound'}, obj.lb = verify_vector(value, 'lb');
                    case {'ub' 'upper' 'upper_bound'}, obj.ub = verify_vector(value, 'ub');
                        
                    % Non-linear constraints
                    case 'nonlcon'   
                        if ~isempty(value)
                            obj.nonlcon{end+1} = verify_function_handle(value,...
                                                                        'Non-linear constraint function');
                        end

                    % Integer constraints
                    case 'intcon'
                        
                        % TODO: intcon can be
                        % - a logical array
                        % - an array of indices
                        % - a cell array with either one of the above, or
                        %   an array of valid integers for that variable                        
                        obj.intcon = value;


                    % Non-linear constraints are evaluated inside the 
                    % objective function
                    case {'constraints_in_objective_function' 'constrinobj'}
                        
                        value = check_logical(value,...
                                              'constraints_in_objective_function');
                        obj.nonlinear_constraints_from_objective_function = value;
                        
                    % The objective function returns multiple arguments,
                    % which are interpreted as the values of mulitple objective 
                    % functions.
                    case {'number_of_objectives' 'objectives'}
                        
                        obj.number_of_objectives = verify_scalar(value, ...
                                                                 'number_of_objectives');                        
                        obj.defines_multiple_objectives = value > 1;                        

                    % Unsupported parameter
                    otherwise
                        warning([mfilename('class') ':unsupported_parameter'],...
                                'Unsupported parameter: "%s"; ignoring...',...
                                parameter);
                end

            end
            
            % Further checks and initializations
            obj.postConstructionChecks();
            
            
            % Validators -------------------------------------------------------
            
            % Check that given parameter is an array of real and finite numbers
            function a = check_double_array(a, parameter_name)
                assert(isnumeric(a) && all(isfinite(a(:))) && all(isreal(a(:))),...
                       [mfilename('class') ':datatype_error'], ...
                       'Argument "%s" must be an array of real, finite values.',...
                       parameter_name);                
            end
            
            function s = verify_scalar(s, parameter_name)
                check_double_array(s, parameter_name);
                assert(isempty(s) || isscalar(s),...
                       [mfilename('class') ':datadims_error'],...
                       'Argument "%s" must be a scalar.',...
                       parameter_name);
            end
            
            function v = verify_vector(v, parameter_name)
                check_double_array(v, parameter_name);
                assert(isempty(v) || isvector(v),...
                       [mfilename('class') ':datadims_error'],...
                       'Argument "%s" must be a vector.',...
                       parameter_name);
            end
            
            function M = verify_matrix(M, parameter_name)
                check_double_array(M, parameter_name);
                assert(isempty(M) || ismatrix(M),...% > R2010b
                       [mfilename('class') ':datadims_error'],...
                       'Argument "%s" must be a matrix.',...
                       parameter_name);
            end
            
            function F = verify_function_handle(F, parameter_name)
                assert(isa(F, 'function_handle'),...
                       [mfilename('class') ':datatype_error'],...
                       'Argument "%s" must be specified with a function handle.',...
                       parameter_name);
            end
            
            % Check that given parameter is a boolean
            function L = check_logical(L, parameter_name)
                if ischar(L)
                    switch lower(L)
                        case {'no' 'off' 'none' 'nope' 'false' 'n'}
                            L = false;
                        case {'yes' 'yup' 'true'  'y'}
                            L = true;
                        otherwise
                            error([mfilename('class') ':datatype_error'], [...
                                  'When specifying argument "%s" via string, ',...
                                  'that string must equal either "yes" or "no".'],...
                                  parameter_name);
                    end
                else
                    assert(islogical(L) && isscalar(value),...
                          [mfilename('class') ':datatype_error'], ...
                          'Argument "%s" must be a logical scalar.',...
                          parameter_name);
                end
            end
            
        end
        
    end
    
    % Operator overloads
    methods       
        
        % Overload subsref, not for indexing, but for evaluation
        function varargout = subsref(obj, X)
            if numel(X)==1 && strcmp(X.type, '()')
                varargout{1} = obj.evaluate(X.subs{:});
            else
                try
                    [varargout{1:nargout}] = builtin('subsref', obj, X);
                catch ME
                    throwAsCaller(ME);
                end
            end
        end

        % Overload feval, for old algorithms still using that
        function varargout = feval(obj, X)
            varargout{1} = obj.evaluate(X);
        end
        
    end
    
    % Public functionality
    methods
        
        % Public accessors; less technical names versions of the private
        % functions below
        function X = getRealX(obj, X_T)
            X = obj.untransformX(X_T);
        end
        
        function X_T = getFakeX(obj, X)
            X_T = obj.transformX(X);
        end
        
    end
    
    % Manual evaluation (for debugging purposes only, hence "hidden"
    methods (Hidden)
        
        % Direct evaluation of function; for debugging purposes        
        function Y = evaluateDirectly(obj, X)
            
            %assert() % numeric, etc. 
            %assert() % A and X have size mismatch, etc.
            
            % Determine constraint violations
            ineq = all(all( obj.A*X <= repmat(obj.b, 1, size(X,2)) ));
            eqc  = all(all( abs(obj.Aeq*X - repmat(obj.beq, 1, size(X,2))) <= eps ));
            ubnd = all(X <= obj.ub);
            lbnd = all(X >= obj.lb);
            intc = all( round(X(obj.intcon))==X(obj.intcon) );
            
            % Get the value of the objective function(s) and non-linear
            % constraints
            if obj.nonlinear_constraints_from_objective_function
                [Y{1:obj.number_of_objectives}, nlc, nlceq] = obj.funfcn(X);
            else
                [nlc, nlceq] = cellfun(@(F)F(X), obj.nonlcon);
                [Y{1:obj.number_of_objectives}] = obj.funfcn(X);
            end
            
            obj.objfcn_evaluations = obj.objfcn_evaluations + 1;
            obj.has_been_evaluated = true;
            
            nlc   = all(nlc <= 0);
            nlceq = all(abs(nlceq) <= eps);
            
            % Display
            cond = {'false', 'true'};
            feas = {'infeasible', 'feasible'};
            
            fprintf(1, [...
                    'Function value: %f\n\n',...
                    'A·x <= b     : %s\n',...
                    'Aeq·x == beq : %s\n',...
                    'x <= UB      : %s\n',...
                    'x >= LB      : %s\n',...
                    'C(X) <= 0    : %s\n',...
                    'Ceq(X) == 0  : %s\n',...
                    'X(intcon) C Z: %s\n',...
                    'Conclusion   : %s\n\n'],...
                    Y{:},...
                    cond{ineq+1},...
                    cond{eqc+1},...
                    cond{ubnd+1},...
                    cond{lbnd+1},...
                    cond{nlc+1},...
                    cond{nlceq+1},...
                    cond{intc+1},...
                    feas{ (ineq && eqc && ubnd && lbnd && nlc && nlceq && intc) + 1});
            
        end
               
    end

    
    % Methods for internal use
    methods (Access = private)
        
        % Validations 
        % ======================================================================
        
        % Elaborate error trapping
        function postConstructionChecks(obj)
            
            assert(~isempty(obj.funfcn),...
                   [mfilename('class') ':function_not_defined'],...
                   'GODLIKE requires at least one objective function.');
               
               
            % Check & adjust variables
            % ------------------------------------------------------------------
                        
            % Ensure tautologies
            obj.is_constrained = ~isempty(obj.A)       || ~isempty(obj.Aeq) || ...
                                 ~isempty(obj.lb)      || ~isempty(obj.ub)  || ...
                                 ~isempty(obj.nonlcon) || ~isempty(obj.intcon);
                             
            if obj.is_constrained
                
                % Checks
                assert(~xor(isempty(obj.A),isempty(obj.b)), [...
                       [mfilename('class') ':inconsistent_inequality'],...
                       'When specifying linear inequalities, both A and b should be ',...
                       'non-empty.']);
                assert(~xor(isempty(obj.Aeq),isempty(obj.beq)), [...
                       [mfilename('class') ':inconsistent_equality'],...
                       'When specifying linear equalities, both Aeq and beq should be ',...
                       'non-empty.']);

                % Clean up the linear constraints                                  
                [obj.A  , obj.b  ] = prepare_linear_constraints(obj.A  , obj.b  , 'inequality');
                [obj.Aeq, obj.beq] = prepare_linear_constraints(obj.Aeq, obj.beq, 'equality'  );
                

                % Clean up the bound constraints                
                [obj.lb, obj.ub] = prepare_bound_constraints(obj.lb, obj.ub);
                
                % Cross-check consistency between LB/UB and A/Aeq
                [obj.A  , obj.b  ] = check_bound_vs_linear_constraints(obj.A  , obj.b  , 'inequality');
                [obj.Aeq, obj.beq] = check_bound_vs_linear_constraints(obj.Aeq, obj.beq, 'equality'  );
                                
                % Check if sizes are all consistent acrcoss all constraints
                check_constraint_dimensions(obj.lb, obj.A, obj.Aeq);
                        
                % Prepare some handydandy indexing variables
                if ~isempty(obj.lb) || ~isempty(obj.ub)
                    
                    obj.nf_lb   = ~isfinite(obj.lb);  
                    obj.nf_ub   = ~isfinite(obj.ub); 

                    obj.fix_var =  obj.lb == obj.ub;                     

                    obj.lb_only = ~obj.nf_lb &  obj.nf_ub & ~obj.fix_var;
                    obj.lb_ub   = ~obj.nf_lb & ~obj.nf_ub & ~obj.fix_var; 
                    obj.ub_only =  obj.nf_lb & ~obj.nf_ub & ~obj.fix_var;  

                    obj.unconst =  obj.nf_lb &  obj.nf_ub & ~obj.fix_var;

                end

            end
            
            % Small helper functions to make life easier
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            
            % Check & clean up linear constraints
            function [M,v] = prepare_linear_constraints(M,v,...
                                                        constraint_type)
                
                if ~isempty(M)
                    
                    switch constraint_type
                        case 'inequality', M_str = 'A';   v_str = 'b'; 
                        case 'equality'  , M_str = 'Aeq'; v_str = 'beq'; 
                    end
                    
                    % Basic assertions
                    assert(isvector(v) && size(M,1)==numel(v),...
                           [mfilename('class') ':invalid_linear_constraints'], [...
                           '%s constraints are implemented as %s·x < %s, where %s and %s ',...
                           'are real and numeric arrays with compatible sizes.'],...
                           constraint_type,...
                           M_str, v_str,...
                           M_str, v_str);
                       
                    % Render into column vector for predictable size
                    v = v(:);
                    
                    % Remove all-zero rows
                    zero_rows = all(M==0,2);
                    if any(zero_rows)

                        assert(all(v(zero_rows) >= 0),...
                               [mfilename('class') ':inconsistent_linear_inequality'], [...
                               'Some rows in the linear %s constraint matrix %s ',...
                               'are all-zero, while the corresponding entries in %s are ',...
                               'non-zero and negative; this constraint is impossible ',...
                               'to meet.'],...
                               constraint_type,...
                               M_str,...
                               v_str);

                        M(zero_rows,:) = [];
                        v(zero_rows  ) = [];
                    end

                    % Remove non-unique rows from A, Aeq 
                    [M_sorted, inds] = sortrows(M);
                    M_repeated_rows  = all(diff(M_sorted)==0, 2);
                    M_repeated_rows  = [false; M_repeated_rows] | [M_repeated_rows; false];
                    
                    if any(M_repeated_rows)
                        
                        v_sorted = v(inds);                        
                        v_repeated = diff( v_sorted( M_repeated_rows ) );
                        assert(all(abs(v_repeated(1:2:end)) <= eps),...
                               [mfilename('class') ':inconsistent_linear_constraint'], [...
                               'Linear %s constraint matrix %s has repeated rows, ',...
                               'but the corresponding elements in %s differ; this ',...
                               'constraint is impossible to meet.'],...
                               constraint_type,...
                               M_str,...
                               v_str);
                        
                        inds(M_repeated_rows)     = [];
                        M_sorted(M_repeated_rows) = [];     M = M_sorted(inds);
                        v_sorted(M_repeated_rows) = [];     v = v_sorted(inds);
                               
                    end
                end
                
            end            
            
            % Check & clean up bound constraints
            function [lb,ub] = prepare_bound_constraints(lb,ub)
                
                % Empties
                if isempty(lb) &&  isempty(ub), return; end                
                if isempty(lb) && ~isempty(ub), lb = -inf(size(ub)); end
                if isempty(ub) && ~isempty(lb), ub = +inf(size(lb)); end
                
                % One bound is scalar, the other is not: replicate
                if isscalar(lb) && ~isscalar(ub)
                    lb = lb*ones(size(ub)); end
                
                if isscalar(ub) && ~isscalar(lb)
                    ub = ub*ones(size(lb)); end
                
                % Basic assertions
                assert(numel(lb) == numel(ub),...
                       [mfilename('class') ':bounds_dimension_mismatch'], ...
                       'Upper and lower bounds differ in size.');
                
                assert(all(lb(:) <= ub(:)),...
                       [mfilename('class') ':lb_exceeds_ub'], [...
                       'All entries in [lb] must be smaller than the corresponding ',...
                       'entries in [ub].']);
                   
                assert(~all(ub(:)==lb(:)),...
                       [mfilename('class') ':bounds_are_equal'], ...
                       'Upper and lower bounds are equal; nothing to do.');
                   
                % Save original size, but change the bound constraints into
                % column vectors
                original_lb_size = size(obj.lb);   obj.lb = obj.lb(:);
                original_ub_size = size(obj.ub);   obj.ub = obj.ub(:);

                % LB/UB may be used to specify dimensions of the decision variable
                if ~isequal(original_lb_size, original_ub_size)
                    warning([mfilename('class') ':lb_ub_dimension_mismatch'],...
                            'Sizes of LB and UB arrays differ; using vector...');
                    original_lb_size = size(obj.lb);
                end
                
                % Set the X-sizes if they havent been specified via option 
                if isempty(obj.original_X_size)
                    obj.original_X_size    = original_lb_size;                
                    obj.transformed_X_size = size(obj.transformX(obj.lb));
                else
                    % TODO: (Rody Oldenhuis) check!
                end
                   
            end
            
            % Check dimensions across all given constraints
            function check_constraint_dimensions(lb, A, Aeq)
                
                have_bounds       = ~isempty(lb);
                have_inequalities = ~isempty(A);
                have_equalities   = ~isempty(Aeq);
                
                % Quick exit
                if sum([have_bounds have_inequalities have_equalities])==1
                    return; end
                
                % Check whether all dimensions of all constraints agree
                if have_bounds 
                    
                    if have_inequalities
                        assert(numel(lb)==size(A,2),...
                               [mfilename('class') ':dimension_mismatch'], [...
                               'Number of elements in bound constraints LB/UB ',...
                               'does not match the number of columns in the linear ',...
                               'inequality constraint matrix A.']);
                    end
                    if have_equalities
                        assert(numel(lb)==size(Aeq,2),...
                               [mfilename('class') ':dimension_mismatch'], [...
                               'Number of elements in bound constraints LB/UB ',...
                               'does not match the number of columns in the linear ',...
                               'equality constraint matrix A_eq.']);
                    end
                
                else % means: we have lb=[], A=nonempty, Aeq=nonempty
                    
                    assert(numel(obj.lb)==size(obj.Aeq,2),...
                           [mfilename('class') ':dimension_mismatch'], [...
                           'Dimension mismatch between linear equality constraint ',...
                           'matrix A_eq and linear inequality constraint matrix ',...
                           'A.']);
                       
                    % Set the X-sizes if they havent been specified yet
                    if isempty(obj.original_X_size)
                        obj.original_X_size = [size(obj.A,2) 1];
                        obj.transformed_X_size = size(obj.transformX(obj.lb));
                    else
                        % TODO: (Rody Oldenhuis) check!
                    end
                    
                end
                
            end
            
            % Check LB/UB and A/Aeq consistency
            function [M,v] = check_bound_vs_linear_constraints(M,v,...
                                                               constraint_type)
                
                if ~isempty(M)
                    
                    % Detect bound constraints hidden in the linear 
                    % (in)equalities. Move them to the bounds if they are
                    % consistent. 
                    
                    warning_issued = false;
                    removals       = false(size(M,1),1);
                    for ii = 1:size(M,1)

                        Arow = M(ii,:);
                        brow = v(ii);
                        
                        % Detect all-zero-but-one rows
                        nz_elements = find(Arow);
                        if numel(nz_elements)==1 
                            
                            if ~warning_issued
                                warning_issued = true;
                                warning([mfilename('class') ':linear_constraint_as_bound'], [...
                                        'At least one linear %s constraint is used as bound ',...
                                        'constraint; use bound constraints for improved ',...
                                        'reliability and efficiency. Re-defining as ',...
                                        'bound constraint...'],...
                                        constraint_type);                                    
                            end
                           
                            % Use the scaled value for comparisons
                            brow = brow/Arow(nz_elements);

                            % Bounds may not be defined yet at this point 
                            if isempty(obj.lb)
                                obj.lb = -inf(obj.original_X_size); end
                            if isempty(obj.ub) 
                                obj.ub = +inf(obj.original_X_size); end

                            % Remove the constraint from the linear (in)
                            % equality constraint and into the bound constraints
                            switch lower(constraint_type)
                                case 'inequality'
                                    obj.ub(ii) = min(obj.ub(ii), brow);
                                    
                                case 'equality'                                     
                                    assert(brow < obj.ub(ii) && brow > obj.lb(ii),...
                                           [mfilename('class') ':inconsistent_constraints'], [...
                                           'What''s worse; the equality constraint mandates ',...
                                           'a value outside the LB/UB domain; inconsistent ',...
                                           'constraints can never be met.']);
                                       
                                    obj.lb(ii) = brow;
                                    obj.ub(ii) = brow;
                            end
                            
                            removals(ii) = true;
                            
                        end
                    end
                    
                    % (The actual deletion)
                    if any(removals)
                        
                        M(removals,:) = [];
                        v(removals)   = [];
                        
                        % (Prevent "0x1 matrix" display)
                        if numel(M)==0
                            M = []; v = []; end
                        
                    end
                    
                end
            end   
            
        end 
        
        
        % Transformations 
        % ======================================================================
        
        % Evaluate transformed function & apply penalties
        function Y = evaluate(obj, X_T)
            
            X = obj.untransformX(X_T);
            Y = obj.funfcnP(X);
            
            % Bookkeeping
            obj.true_X_at_last_evaluation        = X;
            obj.transformed_X_at_last_evaluation = X_T;
            obj.transformed_Y_at_last_evaluation = Y;  
            obj.has_been_evaluated               = true;
            
        end
        
        
        % Transformations 
        % ======================================================================
        
        % Transform possibly-constrained decision variable to its 
        % unconstrained counterpart
        function X = transformX(obj, Z, varargin)
            
            X = Z(:);                                                              
            
            X(obj.lb_only) = sqrt(     Z(obj.lb_only) - obj.lb(obj.lb_only));      
            X(obj.ub_only) = sqrt(obj.ub(obj.ub_only) -      Z(obj.ub_only));      
            X(obj.lb_ub)   = real( asin( 2*(Z(obj.lb_ub) - obj.lb(obj.lb_ub))./ ...
                                   (obj.ub(obj.lb_ub) - obj.lb(obj.lb_ub)) - 1) ); 
                               
            % Fix any unconstrained variables   
            X(obj.fix_var) = [];
            
        end
        
        % Transform unconstrained decision variable to its
        % possibly-constrained counterpart
        function Z = untransformX(obj, X, varargin)
            
            % NOTE: if the bounds are not set, the size of X will not change
            if isempty(obj.original_X_size)
               obj.original_X_size    = size(X); 
               obj.transformed_X_size = size(X); 
            end
               
            Z = zeros(obj.original_X_size);
            Z = Z(:); 
            
            if ~isempty(obj.lb)
                            
                % First insert fixed values...            
                Z( obj.fix_var) = obj.lb(obj.fix_var);
                Z(~obj.fix_var) = X;

                % Handle the integer constraints
                Z(obj.intcon) = round(Z(obj.intcon));

                % Handle the bound constraints
                Z(obj.lb_only) = obj.lb(obj.lb_only) + X(obj.lb_only).^2;
                Z(obj.ub_only) = obj.ub(obj.ub_only) - X(obj.ub_only).^2;
                Z(obj.fix_var) = obj.lb(obj.fix_var);
                Z(obj.unconst) = X(obj.unconst);
                Z(obj.lb_ub  ) = obj.lb(obj.lb_ub) + (obj.ub(obj.lb_ub)-obj.lb(obj.lb_ub)) .* ...
                                                     (sin(X(obj.lb_ub)) + 1)/2; 
                                      
            else
                % Only handle the integer constraints in this case      
                assert(isequal(prod(obj.original_X_size), numel(X)),...
                       [mfilename('class') ':size_mismatch'],...
                       'Expected array with %d elements; received %d elements.',...
                       prod(obj.original_X_size),...
                       numel(X));
                   
                Z(:) = X;
                Z(obj.intcon) = round(Z(obj.intcon));
            end
                                      
            % Reshape
            Z = reshape(Z, obj.original_X_size);
            
        end
        
        % Penalize function according to constraint violation
        function P_fval = funfcnP(obj, x)
                                  
            c   = [];  
            ceq = []; 

            % Evaluate objective function(s)
            if ~obj.nonlinear_constraints_from_objective_function
                [obj_fval{1:obj.number_of_objectives}] = obj.funfcn(x); 
                
            else                    
                arg_out      = cell(1, obj.nonlinear_constraints_from_objective_function+1);
                [arg_out{:}] = obj.funfcn(x);                  
                obj_fval     = arg_out(1:obj.number_of_objectives);
                c            = arg_out{nonlconFcn_in_objFcn + 0};
                ceq          = arg_out{nonlconFcn_in_objFcn + 1};
                
                obj.conffcn_evaluations = obj.conffcn_evaluations + 1;

            end
            
            % Multiple objectives
            obj_fval = [obj_fval{:}];
            
            % Column vector needed for the remainder
            x = x(:);
            
            % Bookkeeping
            obj.true_Y_at_last_evaluation = obj_fval;            
            obj.objfcn_evaluations = obj.objfcn_evaluations + 1;              
                        
            % Initially, we are optimistic
            linear_eq_Penalty   = 0;    
            linear_ineq_Penalty = 0;    
            nonlin_eq_Penalty   = 0;    
            nonlin_ineq_Penalty = 0; 
            integer_Penalty     = 0;

            % Penalize the linear equality constraint violation 
            % required: Aeq*x = beq   
            if ~isempty(obj.Aeq)
                lin_eq    = abs( obj.Aeq*x - obj.beq );                
                sumlin_eq = sum(lin_eq( lin_eq > obj.tolcon ));
                linear_eq_Penalty = penalize_value(sumlin_eq);
            end

            % Penalize the linear inequality constraint violation 
            % required: A*x <= b
            if ~isempty(obj.A)
                lin_ineq = obj.A*x - obj.b;                     
                lin_ineq(lin_ineq <= obj.tolcon) = 0;
                linear_ineq_Penalty = penalize_value(sum(lin_ineq(:)));
            end

            % Penalize the non-linear constraint violations
            % required: ceq = 0 and c <= 0
            if ~isempty(obj.nonlcon)
                
                for ii = 1:numel(obj.nonlcon)
                    
                    [c_, ceq_] = obj.nonlcon{ii}(x); 
                    
                    c   = [c;   c_(:)  ];
                    ceq = [ceq; ceq_(:)];
                end
                
                obj.confcn_evaluations = obj.confcn_evaluations + 1; 
                
            end

            % Process non-linear inequality constraints
            if ~isempty(c)
                violated_c = c > obj.tolcon;                
                nonlin_ineq_Penalty = penalize_value( sum(c(violated_c)) );
            end

            % Process non-linear equality constraints
            if ~isempty(ceq)
                ceq = abs(ceq);
                violated_ceq = (ceq >= obj.tolcon);                 
                nonlin_eq_Penalty = penalize_value( sum(ceq(violated_ceq)) );
            end

            % Return penalized function value
            P_fval = bsxfun(@plus, obj_fval, ...
                     linear_eq_Penalty + linear_ineq_Penalty + ...
                     integer_Penalty + ...
                     nonlin_eq_Penalty + nonlin_ineq_Penalty); 
                 
                 
            % Compute deserved penalties
            function fP = penalize_value(violation)
                
                % Initialize
                persistent exp50
                if isempty(exp50)
                    exp50 = exp(50); end

                % Safety condition
                if violation <= obj.tolcon
                    fP = 0; 
                    return; 
                end

                % Scaling parameter    
                scale = min(1e60, violation/obj.tolcon);          
                
                % Linear penalty to avoid overflow
                if scale*violation > 50
                    fP = exp50*(1 + scale*violation) - 1;

                % Exponential penalty otherwise
                else
                    fP = exp(scale*violation) - 1;
                end

            end % Penalize


        end 

    end

end




