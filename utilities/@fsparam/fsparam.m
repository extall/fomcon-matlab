classdef fsparam
    %FSPARAM Fractional-order system simulation parameters
    %
    % Creates a new simulation parameters structure
    % used for FOTF system simulation using one of the following
    % approximations: Oustaloup filter, refined Oustaloup filter
    %
    % Usage: P = FSPARAM(PLANT, APPROX, W, N)
    % 
    % where
    %           PLANT - FOTF object (plant) OR
    %                   Symbolic expression in s which evaluates to
    %                   a valid (parametrized) transfer function
    %           APPROX - approximation method, can be
    %                    either 'gl', 'oust', or 'ref'
    %           W - approximation frequency range in form [WB; WH] (rad/s)
    %           N - approximation order
    %
    % Any parameter except G (FOTF object) can be omitted in which case
    % default values are used.
    %
    % Note that the given FOTF object is stored as a copy, not a reference,
    % i.e. you will need to update the structure manually when the initial 
    % object is changed.
    %
    % Defaults: P = FSPARAM(..., 'oust', [1e-3;1e3], 5);
    %
    % See also: oustapp, fid, fpid_optimize, 
    
    properties
        plant           % FOTF object
        fotf_expr       % Symbolic expression corresponding to the plant
        approx          % Approximation method
        w               % Approximation valid frequency range [rad/s]
        N               % Approximation order
    end
    
    methods
        
        % Initializer
        function p = fsparam(G, method, w, N)
           
            % Defaults
            d_method = 'oust';
            d_w = [0.001, 1000];
            d_N = 5;
            
            % Check input argument number
            if nargin == 0
                error('Not enough input arguments.');
            end
            
            if nargin < 4
               N = d_N;
            end
            
            if nargin < 3
               w = d_w;
            end
                     
            if nargin < 2
               method = d_method;
            end
            
            % Check input arguments
            
            % Method
            method = lower(method);
            switch(method)
                case {'gl', 'oust', 'ref'}
                    % Do nothing
                
                otherwise
                    method = 'oust';
                    warning('FSPARAM:BadApproximationMethod', ...
                    'Invalid approximation method supplied, using ''oust''');
            end
            
            % Frequency
            if w(1) > w(2)
                
                temp = w(1);
                w(1) = w(2);
                w(2) = temp;
                
                warning('FSPARAM:FrequenciesSwapped', ...
                'Frequency bounds are such that wb > wh; swapping');
                
            end
            
            % Check approximation order
            if N < 1
               
                N = 1;
                warning('FSPARAM:BadApproximationOrder', ...
                'Approximation order must be N>=1, setting to 1');
               
            end
            
            % Set parameters
            
            % Plant/expression
            if (isa(G, 'fotf'))
                p.plant = G;
                p.fotf_expr = '';
            elseif (isa(G, 'char'))
                p.plant = [];
                p.fotf_expr = G;
            else
                error('Wrong type of system or expression given!');
            end
            
            p.approx = lower(method);
            p.w = w;
            p.N = N;
            
        end
        
    end
    
end

