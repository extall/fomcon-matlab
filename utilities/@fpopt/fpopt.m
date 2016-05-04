classdef fpopt
    %FPOPT Creates a new fractional PID optimization parameter structure
    %
    % Usage: OPT = FPOPT(TYPE,P,PMAX,PMIN,METRIC,ALG,SP,MARGINS,SENS,ULIM,
	%                    WGT,GAINVAR,SIMTIME,CANCELZERO,STRICT,OPTOP)
    %
    % where:    
    %            TYPE       - optimization type, can be set to fix 'n' (none), 
    %                         'e' (exponents) or 'g' (gains) for fixing
    %                         exponents or gains
    %            P          - parameters to optimize: [Kp Ki Kd lam mu]
    %            PMAX       - maximum allowed values of P
    %            PMIN       - minimum allowed values of P
    %            METRIC     - fractional control system performance metric,
    %                         can be 'ise', 'iae', 'itse' or 'itae'
	%			 ALG        - optimization algorithm, can be 'nelder-mead',
	%						  'interior-point', 'sqp' or 'active-set' - the
	%                         last three choices require the availability of
	%                         the Optimization toolbox
    %            SP         - set-point for optimization
    %            MARGINS    - gain and phase margin specifications in
    %                         form [GM, gmExact; PM, pmExact], GM in [dB],
    %                         PM in deg and gmExact and pmExact are boolean
    %                         values requiring exact specification to be
    %                         met (if set to true); otherwise optimization
    %                         will consider these parameters as minimum
    %                         requirements
    %            SENS    	- sensitivity function specifications in form
    %                         [wt, A;   % frequency wt [rad/s] and A [dB]
    %                          ws, B]   % frequency ws [rad/s] and B [dB]
    %            ULIM    	- control signal limits in form [UMIN; UMAX]
    %            WGT        - if nonempty, will be taken as the control law
    %                         weight in the optimization objective
    %                         function. If empty, the objective function
    %                         remains unchanged.
    %            GAINVAR    - robustness to gain variation settings in form
    %                         [WCG; WRANGE] (rad/s). The following options
    %                         are additionally possible:
    %                         * Set WCG=WRANGE, and the optimizer will look
    %                           for a local extremal of dPhase/dw at WCG.
    %                         * Set WRANGE = -1 and ONLY the WCG
    %                           specification will be considered. This is
    %                           useful when the bandwidth of the system has
    %                           to be set precisely.
    %            SIMTIME    - simulation time step min/max and max
    %                         simulation time (seconds) in form [DTMIN;
    %                         DTMAX; MAXTIME]
	%			 CANCELZERO - (true/false) enable zero cancelation for non-proper systems
    %            STRICT     - (true/false) strict specification requirement (including i.c.)
    %            OPTOP      - special options for optimization algorithm
    %
    % Note: call without arguments for default values. You may use []
    %       instead of parameters MARGIN, SENS and ULIM to optimize without regard
	%		to these specifications.
    %
    % See also: fpid_optimize, margin, csens
    
    properties
        type        % Optimization type
        
        p           % Optimized parameters
        pmax        % Max p values
        pmin        % Min p values
        
        metric      % Performance metric
		alg 		% Optimization algorithm
        
        sp          % Set-point
        
        margins     % Gain and phase margin specifications
        sens        % Sensitivity & compl. sensitivity function parameters
        ulim        % Control signal limits
        wgt         % Control signal weight in objective function
        
        gainvar     % Robustness to gain variation parameters
        
        simtime     % Maximum simuation time
		
		cancelzero	% Zero cancelation for non-proper systems

        strict      % Strict option, including initial conditions
        
        optop     	% Optimization options
        
    end
    
    methods
        
        function opt = fpopt(type, p, pmax, pmin, metric, alg, sp, ...
                             margins, sens, ulim, wgt, gainvar, simtime, ...
                             cancelzero, strict, optop)
            
            % Default values
            type_d       = 'n';
            p_d          = ones(1,5);
            pmax_d       = [1000 1000 1000 2 2];
            pmin_d       = [-1000 -1000 -1000 0.01 0.01];
            metric_d     = 'ise';
			alg_d        = 'nelder-mead';
            sp_d         = 1;
            mar_d        = [10, 0; 60, 0];
            sens_d       = [100, -20; 0.01, -20];
            ulim_d       = [-500; 500];
            wgt_d        = [];
            gainvar_d    = [0.01; 0.005];
            simtime_d    = [0.001; 0.5; 10000];
			cancelzero_d = false;
            strict_d     = false;
            optop_d      = optimset('MaxIter', 50);
            
            % Check number of input arguments
            if nargin < 16
                optop = optop_d; 
            end
            
            if nargin < 15
                strict = strict_d;
            end
			
            if nargin < 14
				cancelzero = cancelzero_d;
            end
            
            if nargin < 13
                simtime = simtime_d;
            end
            
            if nargin < 12
                gainvar = gainvar_d;
            end
            
            if nargin < 11
                ulim = ulim_d;
            end
            
            if nargin < 10
                wgt = wgt_d;
            end
            
            if nargin < 9
                sens = sens_d;
            end
            
            if nargin < 8
                margins = mar_d;
            end
            
            if nargin < 7
                sp = sp_d;
            end
			
			if nargin < 6
				alg = alg_d;
			end
            
            if nargin < 5
                metric = metric_d;
            end
            
            if nargin < 4 || isempty(pmin)
               pmin = pmin_d; 
            end
            
            if nargin < 3 || isempty(pmax)
               pmax = pmax_d; 
            end
            
            if nargin < 2 || isempty(p)
               p = p_d; 
            end
            
            if nargin == 0
               type = type_d; 
            end
            
            % Check parameters
            type = lower(type);
            switch(type)
                case 'n'
                case 'e'
                case 'g'
                    % Do nothing
                otherwise
                    type = type_d;
                    warning('FPOPT:BadType', 'Invalid optimization type, using default');
            end
            
            metric = lower(metric);
            switch(metric)
               
                case 'ise'
                case 'iae'
                case 'itse'
                case 'itae'
                    % Do nothing
                otherwise
                    metric = metric_d;
                    warning('FPOPT:BadMetric', 'Invalid metric type, using default');
            end
            
            % Constrain the weight to the interval (0, 1).
            if ~isempty(wgt)
                if (wgt < 0) || (wgt > 1)
                    error('FPOPT:ControlLawWeightOutsideAcceptableRange', ...
                          'Control law weight is outside the acceptable range');
                end
            end
            
            % Set parameters
            opt.type       = type;
            opt.p          = p;
            opt.pmax       = pmax;
            opt.pmin       = pmin;
            opt.metric     = metric;
			opt.alg        = alg;
            opt.sp         = sp;
            opt.margins    = margins;
            opt.sens       = sens;
            opt.ulim       = ulim;
            opt.wgt        = wgt;
            opt.gainvar    = gainvar;
            opt.simtime    = simtime;
			opt.cancelzero = cancelzero;
            opt.strict     = strict;
            opt.optop      = optop;
            
        end
        
    end
    
end

