classdef fpid_simopt
    %FPID_SIMOPT Simulink simulation options structure for FPID_OPTIMIZE
    
    properties
        
        modelname       % Model to use for simulation
        
        tmax            % Simulation stop time
        dtmin           % Minimum time step
        dtmax           % Maximum time step
        
        ulim            % Control signal limits
        
        r               % Reference value for control
        
        w               % Approximation frequency range [rad/s]
        N               % Appriximation order
        type            % Approximation type
        
        cancelzero      % Whether to ensure a proper system
        
    end
    
    methods
        % Initializer
        function s = fpid_simopt(modelname, tmax, dtmin, dtmax, ulim, r, w, N, type, cancelzero)

            % Default values
            modelname_d  = 'fpid_optimize_model.mdl';
            
            tmax_d       = 100;            % Time vector and min/max steps
            dtmin_d      = 0.01;
            dtmax_d      = 0.5;
            
            ulim_d       = [];             % Control signal limits
            
            r_d          = 1;              % Reference value
            
            w_d          = [0.001; 1000];  % Approximation parameters
            N_d          = 5;
            type_d       = 'oust';
            
            cancelzero_d = false;          % Zero cancelation
            
            % Check input argument number, use default values
            if nargin < 10 || isempty(cancelzero), cancelzero = cancelzero_d; end
            if nargin < 9  || isempty(type),       type  = type_d;            end
            if nargin < 8  || isempty(N),          N     = N_d;               end
            if nargin < 7  || isempty(w),          w     = w_d;               end
            if nargin < 6  || isempty(r),          r     = r_d;               end
            if nargin < 5  || isempty(ulim),       ulim  = ulim_d;            end
            if nargin < 4  || isempty(dtmax),      dtmax = dtmax_d;           end
            if nargin < 3  || isempty(dtmin),      dtmin = dtmin_d;           end
            if nargin < 2  || isempty(tmax),       tmax  = tmax_d;            end
            if nargin < 1  || isempty(modelname),  modelname = modelname_d;   end
            
            % Check arguments
            
            % Check model existence, revert to default if file not found
            if strcmpi(modelname, 'default') || isempty(which(modelname))
                modelname = modelname_d; 
            end
            
            
            if dtmin > dtmax                % Time step min/max are swapped
                
                temp = dtmax;
                dtmax = dtmin;
                dtmin = temp;
                
                warning('FPIDSIMOPT:TimeStepsSwapped', ...
                        'Time steps are such that dtmin > dtmax; swapping');
                    
            end
            
            type = lower(type);             % Check approximation type
            switch(type)
                case {'oust', 'ref'}
                    % Do nothing
                otherwise
                    type = 'oust';
                    warning('FPIDSIMOPT:BadApproximationMethod', ...
                            'Invalid approximation method supplied, using ''oust''');
            end
            
            if w(1) > w(2)                  % Frequencies are swapped
                
                temp = w(1);
                w(1) = w(2);
                w(2) = temp;
                
                warning('FPIDSIMOPT:FrequenciesSwapped', ...
                        'Frequency bounds are such that wb > wh; swapping');
            end
            
            if N < 1                        % Check approximation order
               
                N = 1;
                warning('FPIDSIMOPT:BadApproximationOrder', ...
                        'Approximation order must be N>=1, setting to 1');
               
            end
            
            if ~isempty(ulim) && ulim(1) > ulim(2)  % Control signal limits are swapped
                
                temp = ulim(1);
                ulim(1) = ulim(2);
                ulim(2) = temp;
                
                warning('FPIDSIMOPT:ControlSignalLimitsSwapped', ...
                        'Control signal limits are such that umin > umax; swapping');
                    
            end
            
            % Assign object properties
            s.modelname  = modelname;
            s.tmax       = tmax;
            s.dtmin      = dtmin;
            s.dtmax      = dtmax;
            s.ulim       = ulim;
            s.r          = r;
            s.w          = w;
            s.N          = N;
            s.type       = type;
            s.cancelzero = cancelzero;
             
        end
    end
end

