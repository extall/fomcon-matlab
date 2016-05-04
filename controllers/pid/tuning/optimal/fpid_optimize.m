function [Kp, Ki, Kd, lam, mu, Z] = fpid_optimize(fsim, fopt, G, usesim)
%FPID_OPTIMIZE Obtain optimized settings for fractional PID controller
%
%   Usage: [KP, KI, KD, LAM, MU, Z] = FPID_OPTIMIZE(FSIM, FOPT, G, USESIM)
%
%   where  FSIM   - FOTF simulation structure,
%          FOPT   - FPID optimization structure of type FPOPT,
%          G      - optional: any valid MATLAB SISO LTI model (TF, ZPK, SS),
%                   which is used instead of FSIM.PLANT. In this case correct
%                   FSIM parameters still need to be supplied because they
%                   are used for FPID approximation while FSIM.PLANT may be
%                   an empty FOTF object.
%          USESIM - optional: if provided a string with model name/path, a
%                   Simulink model will be used for simulating the control
%                   system time-domain response. This option requires
%                   Simulink to be installed. By default the option is an
%                   empty string (i.e. Simulink will not be used).
%
%          KP, KI, KD, LAM, MU - obtained controller parameters.
%                   
%          Z      - structure with the following parameters:
%                   .initialParams  - initial PID parameters
%                   .timeSpan       - calculated time span for optimization
%                   .finalIndex     - achieved performance index
%                   .finalSpecs     - achieved performance specification
%                                     deviations from the desired values
%
%   See also: fsparam, fpopt, tf, zpk, ss

    % Load configuration parameters
    config = fomcon('config');

    if nargin < 2
        error('FPIDOPTIMIZE:NotEnoughInputArguments', ...
            'Not enough input arguments.');
    end

    % Check number of input arguments    
    if nargin >= 3 && ~isempty(G)
       
        nonfotf = true;
        nfmodel = G;
        
    else
        
        nonfotf = false;
        nfmodel = [];
        
    end
    
    % Do not use Simulink by default
    if nargin < 4 || isempty(usesim)
        usesim = [];
    end
    
    % Get simulation options
    G       = fsim.plant;
    method  = fsim.approx;
    wb      = fsim.w(1);
    wh      = fsim.w(2);
    N       = fsim.N;
    
    % Get optimization parameters
    optType = fopt.type;
    
    % Parameters
	
	% Initial values
    Kp.val  = fopt.p(1);
	Ki.val  = fopt.p(2);
	Kd.val  = fopt.p(3);
	lam.val = fopt.p(4);
	mu.val  = fopt.p(5);
	
	% Max values
	Kp.max  = fopt.pmax(1);
	Ki.max  = fopt.pmax(2);
	Kd.max  = fopt.pmax(3);
	lam.max = fopt.pmax(4);
	mu.max  = fopt.pmax(5);
	
	% Min values
	Kp.min  = fopt.pmin(1);
	Ki.min  = fopt.pmin(2);
	Kd.min  = fopt.pmin(3);
	lam.min = fopt.pmin(4);
	mu.min  = fopt.pmin(5);
	
	% Performance metric
	metric = fopt.metric;
	
	% Optimization algorithm
	alg    = fopt.alg;
	
	% Performance specifications
	margins    = fopt.margins;
	sens       = fopt.sens;
    ulim       = fopt.ulim;
    wgt        = fopt.wgt;
	
    % Time parameters
    mindt      = fopt.simtime(1);
    maxdt      = fopt.simtime(2);
    maxtime    = fopt.simtime(3);
    
    % Gain variation specifications
    gainvar    = fopt.gainvar;
    
	% Set point
	sp  	   = fopt.sp;
	
	% Zero cancelation
	cancelzero = fopt.cancelzero;
    
    % Strictness
    switch fopt.strict
        case 0
            strictOption = [];
        case 1
            strictOption = 'strict';
    end
	
	% Optimization options
	op = fopt.optop;
    
    % Initial controller
    Gc_init = oustapid(Kp.val, Ki.val, lam.val, Kd.val, mu.val, ...
                       wb, wh, N, method);
    
    % Cancel zeros
    if fopt.cancelzero
        Gc_init = toproper(Gc_init, wh);
    end
    
    % Simulate initial system
    if nonfotf
        G_plant = nfmodel;
    else
		G_init  = oustapp(G, wb, wh, N, method);          
        G_plant = G_init;
    end
    
    % If the plant has delays, convert it to state-space
    if ~fleq(G_plant.ioDelay, 0)
        G_plant = ss(G_plant);
    end
    
    % Obtain initial control system response
    [ig1, t] = step(feedback(G_plant*Gc_init,1));
	
    dt = (t(2)-t(1))/2;          % dt
    
    % Check against simulation time parameters
    if dt > maxdt
        dt = maxdt;
    elseif dt < mindt
        dt = mindt;
    end
    
    if t(end) < maxtime          % Ensure stability, considering maxtime
        t=0:dt:t(end)*2;             
    else
        t=0:dt:maxtime;
    end
    
    % Save time span
    Z.timeSpan = t;
    
    % Set options for objective function and non-linear constraints
    
    % Simulation options
    opt.nonfotf = nonfotf;
    opt.nfmodel = nfmodel;
    
    opt.wb      = wb;
    opt.wh      = wh;
    opt.type    = method;
    opt.N       = N;
    
    % Metric
    opt.optim   = metric;
    
    % Time vector/step
    opt.t  = t;
    opt.dt = dt;
    
	% Specifications
	opt.margins = margins;
	opt.sens    = sens;
	opt.ulim    = ulim;
    opt.wgt     = wgt;
    opt.gainvar = gainvar;
    
    % Time parameters
    opt.tmax = maxtime;
    opt.dtmin = mindt;
    opt.dtmax = maxdt;
	
	% Set point
	opt.sp = sp;
	
	% Zero cancelaction
	opt.cancelzero = cancelzero;
    
    % Use Simulink?
    opt.usesim = usesim;
    
    % Check reference signal type
    if (isa(sp,'struct') && isempty(usesim))
        error('FPIDOPTIMIZE:CannotUseCustomRef', ...
             'Unable to use custom reference signal unless Simulink is used');
    end
	
    % Parameters to be sent
    opt.kp  = Kp.val;
    opt.ki  = Ki.val;
    opt.kd  = Kd.val;
    opt.lam = lam.val;
    opt.mu  = mu.val;
    
    % Save initial parameters
    Z.initialParams = [Kp.val Ki.val Kd.val lam.val mu.val];
    
    % Fix parameters
    opt.fix = optType;
    
    % Begin optimization
    
    % Measure elapsed time
    elTime = tic();
    
	% Optimize according to selected algorithm
	switch(lower(alg))
	
		case 'nelder-mead'
	
			% Optimize based on type
			switch(optType)
			   
				case 'n'
					 
					 x_opt = optimize(@(x) fpid_optimfun(x,G,opt), ...
								 [Kp.val Ki.val Kd.val lam.val mu.val], ...
								 [Kp.min Ki.min Kd.min lam.min mu.min], ...
								 [Kp.max Ki.max Kd.max lam.max mu.max], ...
								 [], [], [], [], ...
								 @(x) fpid_optimnlc(x,G,opt), ...
								 strictOption, ...
								 op);
                             
					 % Set gains
					 Kp.val = x_opt(1);
					 Ki.val = x_opt(2);
					 Kd.val = x_opt(3);
					 
					 % Set exponents
					 lam.val = x_opt(4);
					 mu.val = x_opt(5);
						
				case 'e'
					
					x_opt = optimize(@(x) fpid_optimfun(x,G,opt), ...
								[Kp.val Ki.val Kd.val], ...
								[Kp.min Ki.min Kd.min], ...
								[Kp.max Ki.max Kd.max], ...
								[], [], [], [], ...
								@(x) fpid_optimnlc(x,G,opt), ...
								strictOption, ...
								op);
							
					% Set gains
					Kp.val = x_opt(1);
					Ki.val = x_opt(2);
					Kd.val = x_opt(3);
								
				case 'g'

					x_opt = optimize(@(x) fpid_optimfun(x,G,opt), ...
								[lam.val mu.val], ...
								[lam.min mu.min], ...
								[lam.max mu.max], ...
								[], [], [], [], ...
								@(x) fpid_optimnlc(x,G,opt), ...
								strictOption, ...
								op);
							
					% Set exponents
					lam.val = x_opt(1);
					mu.val = x_opt(2);
				
			end
		
		case {'interior-point', 'sqp', 'active-set'}
		
			% Set the internal algorithm
			op.Algorithm = lower(alg);
		
			% Optimize based on type
			switch(optType)
			   
				case 'n'
					 
					 x_opt = fmincon(@(x) fpid_optimfun(x,G,opt), ...
								 [Kp.val Ki.val Kd.val lam.val mu.val], ...
								 [], [], [], [], ...
								 [Kp.min Ki.min Kd.min lam.min mu.min], ...
								 [Kp.max Ki.max Kd.max lam.max mu.max], ...
								 @(x) fpid_optimnlc(x,G,opt), ...
								 op);
                             
					 % Set gains
					 Kp.val = x_opt(1);
					 Ki.val = x_opt(2);
					 Kd.val = x_opt(3);
					 
					 % Set exponents
					 lam.val = x_opt(4);
					 mu.val = x_opt(5);
						
				case 'e'
					
					x_opt = fmincon(@(x) fpid_optimfun(x,G,opt), ...
								[Kp.val Ki.val Kd.val], ...
								[], [], [], [], ...
								[Kp.min Ki.min Kd.min], ...
								[Kp.max Ki.max Kd.max], ...
								@(x) fpid_optimnlc(x,G,opt), ...
								op);
							
					% Set gains
					Kp.val = x_opt(1);
					Ki.val = x_opt(2);
					Kd.val = x_opt(3);
								
				case 'g'

					x_opt = fmincon(@(x) fpid_optimfun(x,G,opt), ...
								[lam.val mu.val], ...
								[], [], [], [], ...								
								[lam.min mu.min], ...
								[lam.max mu.max], ...
								@(x) fpid_optimnlc(x,G,opt), ...
								op);
							
					% Set exponents
					lam.val = x_opt(1);
					mu.val = x_opt(2);
				
			end
	end
	
	elTime = toc(elTime);
    
	% Set new parameters
	Kp 		= Kp.val;
	Ki 		= Ki.val;
	Kd 		= Kd.val;
	lam 	= lam.val;
	mu 		= mu.val;
    
	% Optimized controller
	fpid_o = oustapid(Kp, Ki, lam, Kd, mu, wb, wh, N, method);
    
    % Get core parameters
    numPts = config.Core.Frequency_domain_computations.Num_points;
    
    if nonfotf
        
        % Get Bode data in frequencies of interest
        [Mag, Ph, w] = bode(nfmodel*fpid_o, ...
                            logspace(get10exp(opt.wb),...
                            get10exp(opt.wh),numPts));
    else                
        
        % Get Bode data in frequencies of interest
        [Mag, Ph, w] = bode(G*fracpid(Kp, Ki, lam, Kd, mu), ...
                            logspace(get10exp(opt.wb),...
                            get10exp(opt.wh),numPts));
    end
        
    elTimStr = ['Elapsed time: ' stohms(elTime)];
    
    disp(elTimStr);
    disp(char(ones(1,length(elTimStr))*'-'));

    % Get new gain and phase margins
    [Gm_new, Pm_new] = margin(Mag,Ph,w);
    
	numSigDig = config.Core.General.Model_significant_digits;
    Gm_new = num2str(20*log(Gm_new)/log(10), numSigDig);
    Pm_new = num2str(Pm_new, numSigDig);
    
    disp(['Gain margin [dB]: ' Gm_new]);
    disp(['Phase margin [deg]: ' Pm_new]);
    disp('Optimization completed.');
    
end

