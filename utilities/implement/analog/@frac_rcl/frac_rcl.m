classdef frac_rcl
    %FRAC_RCL Implementation of a circuit exhibiting fractance properties 
    % using discrete R/C/L components
    
    properties (SetAccess = protected)
        
        % Fractional-order transfer function, system gain and frequency
        % range (points where the frequency response must be computed)
        model = fotf;
        K     = [];
        w     = [];
        
        % Fractance implementation structure
        structure = '';
        
        % Approximation type
        implementation = '';
        
        % Circuit type
        circuit_type = '';
        
        % Discrete component values
        R = [];
        C = [];
        L = [];
        
        % Additional parameters
        params  = [];
        results = [];
        
    end
    
    methods
        
        % Initializer
        function fract = frac_rcl(G, st, imp, w, params)
            
            % Set model, structure and type
            fract.model          = G;
            fract.structure      = st;
            fract.implementation = imp;
            
            % Determine type of structure
            str_fun = str2func(st);
            str = str_fun(); str = ('RCL'.*str);
            str(fleq(str,0))=[]; str = char(str);
            fract.circuit_type = str;
            
            % Set frequency range
            fract.w = w;
            fract.params = params;
            
            % Add structure type to parameters
            params.model     = G;
            params.structure = st;
            params.w         = w;
            
            % Implement the system
            imp_fun = str2func(imp);
            [R, C, L, K, results ] = imp_fun(params);
            
            % Set the corresponding R,C and L parameters, 
            % block gain and implementation results
            fract.R       = R;
            fract.C       = C;
            fract.L       = L;
            fract.K       = K;
            fract.results = results;

        end

    end
    
end

