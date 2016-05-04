function [c, ceq] = fpid_optimnlc(x, G, opt)
%FPID_OPTIMNLC Non-linear constraints for FPID_OPTIMIZE
    
    % Load configuration parameters
    config = fomcon('config');

    numPts      = config.Core.Frequency_domain_computations.Num_points;
    phaseWeight = config.FPID_Optimizer.Frequency_domain_computations.Phase_comp_weight;
    wcgWeight   = config.FPID_Optimizer.Frequency_domain_computations.W_cg_comp_weight;
    numPhasePts = config.FPID_Optimizer.Frequency_domain_computations.Num_points_phase;
    gmLimit     = 1000; %config.FPID_Optimizer.Frequency_domain_computations.Gain_margin_limit;
    pmLimit     = 1000; %config.FPID_Optimizer.Frequency_domain_computations.Phase_margin_limit;

    % Get controller
    switch lower(opt.fix)
        case 'e'
			Kp = x(1);
			Ki = x(2); lam = opt.lam;
			Kd = x(3); mu  = opt.mu;
        case 'g'
			Kp = opt.kp;
			Ki = opt.ki; lam = x(1);
			Kd = opt.kd; mu  = x(2);
        case 'n'
			Kp = x(1);
			Ki = x(2); lam = x(4);
			Kd = x(3); mu  = x(5);
        otherwise
            error('Wrong switch for opt.fix!');
    end
    
    % Equality constraints
    ceq = [];
    
    % Inequality constraints
    c = [];
	
    Gc = oustapid(Kp, Ki, lam, Kd, mu, opt.wb, opt.wh, opt.N, opt.type);
    
	% Get open-loop control plant transfer function and
	% closed-loop controller output transfer function
    if opt.nonfotf
        open_loop_cs  = opt.nfmodel*Gc;
    else
        open_loop_cs  = oustapp(G, opt.wb, opt.wh, opt.N, opt.type) * Gc; 
    end
    
    % Gain and phase margins
    if ~isempty(opt.margins)
        
        Gm = opt.margins(1,1);
        Pm = opt.margins(2,1);
        
        % Get Bode data in frequencies of interest
        [Mag, Ph, w] = bode(open_loop_cs, ...
            logspace(get10exp(opt.wb),get10exp(opt.wh),numPts));
    
        % Compute open-loop parameters
        [Gm_real, Pm_real] = margin(Mag,Ph,w);
    
        % Get Gm_real in dB
        Gm_real = 20*log(Gm_real)/log(10);
		
		% Check values
		if isinf(Gm_real), Gm_real = sign(Gm_real) * gmLimit; end
		if isinf(Pm_real), Pm_real = sign(Pm_real) * pmLimit; end
    
        % Compute difference
        GmDif = Gm-Gm_real;
        PmDif = Pm-Pm_real;
        
        % Gain/phase margin constraints (based on Exact requirement)
        if opt.margins(1,2)
            ceq = [ceq; GmDif];
        else
            c   = [c; GmDif];
        end
        
        if opt.margins(2,2)
            ceq = [ceq; PmDif];
        else
            c   = [c; PmDif];
        end
        
    end
	
	% Sensitivity function evaluation
    if ~isempty(opt.sens)
	
		% Complementary sensitivity
		wt = opt.sens(1,1);
		A  = opt.sens(1,2);
		
		% Sensitivity function
		ws = opt.sens(2,1);
		B  = opt.sens(2,2);

		[wt_real, ws_real] = csens(open_loop_cs, A, B, ...
            logspace(get10exp(opt.wb),get10exp(opt.wh),numPts));
        
		% Sensitivity function constraints
		c = [c; wt_real - wt; ws_real - ws];
		
    end
    
    % Robustness to gain variation
    if ~isempty(opt.gainvar)
        wcg   = opt.gainvar(1);
        wcg_h = opt.gainvar(2);
        
        % Compute the necessary frequency-domain parameters
        [Mag, Ph, w] = bode(open_loop_cs, ...
                       logspace(get10exp(opt.wb),get10exp(opt.wh),numPts));
        [ig1, ig2, ig3, wcg_real] = margin(Mag,Ph,w);
        
        % Use a single point instead of a range
        if fleq(wcg,wcg_h)
            
            % If the two parameters are equal, 
            % make d arg(CG) / d w vanish at wcg
            w_comp = logspace(get10exp(opt.wb),get10exp(opt.wh),numPhasePts);
            phase  = rad2deg(angle(squeeze(freqresp(open_loop_cs, w_comp))));
            dCG    = diff(phase);
            ind    = find(w_comp>=wcg_real, 1)-1;

            c      = [c; wcgWeight*(wcg-wcg_real)^2; phaseWeight*(abs(dCG(ind)))];
            
        % Consider only the wcg specification
        elseif ~fleq(wcg,0) && fleq(wcg_h,-1)
            c      = [c; wcgWeight*(wcg-wcg_real)^2];
            
        % Consider the range
        else
            % If wcg > wcg_h then swap
            if wcg > wcg_h, temp = wcg; wcg = wcg_h; wcg_h = temp; end
            
            % Obtain range
            wcg_l = wcg^2 / wcg_h;
            
            
            deg_flat = angle(squeeze(freqresp(open_loop_cs, ...
                logspace(log10(wcg_l), ...
                log10(wcg_h), ...
                numPhasePts))));
            
            flatns = sum(diff(rad2deg(deg_flat)).^2);
            
            % Flatness of the phase curve near critical frequency wcg
            c  = [c; wcgWeight*(wcg-wcg_real)^2; phaseWeight*flatns];
        end
        
    end

end

