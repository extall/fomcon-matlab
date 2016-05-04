function z = fpid_optimfun(x,G,opt)
%FPID_OPTIMFUN Optimal FPID Design helper function
%        Usage: z = fpid_optimfun(x,G,opt)
%        where              
%               x(1) = Kp   or  x(1) = lam (gains fixed)
%               x(2) = Ki       x(2) = mu  
%               x(3) = Kd
%               ---------- (exponents fixed)
%               x(4) = lam
%               x(5) = mu
%               
%               G - FOTF object, plant transfer function
%               opt.optim - optimization metric: 'ise', 'iae', 'itse' or
%                           'itae'.
%               opt.t, opt.dt - time vector, time step
%               opt.wb, opt.wh, opt.N, opt.type - Oustaloup filter params,
%                                                 type: 'oust' or 'ref'.
%               opt.fix - 'e' to fix exponents, 'g' to fix gains, 'n' to
%                         fix nothing; supply fixed parameters as:
%                         opt.kp, opt.ki, opt.kd, opt.lam, opt.mu
%               opt.nonfotf - (true/false): is a LTI model of type TF, ZPK,
%                             SS is used instead of a FOTF object?
%               opt.nfmodel - corresponding LTI model
%				opt.sp - set point
%				opt.cancelzero - (true/false) whether to cancel zeros
%				                 to obtain a proper model
%               opt.usesim - string; if unempty, simulate specified model
%               in Simulink, use 'fpid_optimize_model.mdl' for the
%               default model. If empty, Simulink is not invoked.


    % Load configuration parameters
    config = fomcon('config');
    
    indexWeight = config.FPID_Optimizer.Time_domain_computations.Performance_index_weight;

    if nargin < 3
        error('Not enough input arguments.');
    end
    
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
    
    if opt.nonfotf
        Z_c  = oustapid(Kp, Ki, lam, Kd, mu, opt.wb, opt.wh, opt.N, opt.type);
        Z_p  = opt.nfmodel;
    else
        G1   = oustapp(G, opt.wb, opt.wh, opt.N, opt.type);
        Z_c  = oustapid(Kp, Ki, lam, Kd, mu, opt.wb, opt.wh, opt.N, opt.type);
        Z_p  = G1;
    end
    
    % Convert to proper objects?
    if opt.cancelzero
        Z_c = toproper(Z_c, opt.wh);
        Z_p = toproper(Z_p, opt.wh);
    end
    
    % Handle delay systems in case of linear simulation
    if ~isa(Z_p,'ss') && ~isempty(Z_p.ioDelay) && Z_p.ioDelay ~= 0
        Z_p = ss(Z_p);
    end

    % Time vector
    t  = opt.t;
    dt = opt.dt;
	
    % Determine simulation type
    if ~isempty(opt.usesim)
        % Use Simulink
        
        % Generate simulation options structure
        sopt = fpid_simopt(opt.usesim, opt.tmax, opt.dtmin, opt.dtmax, ...
                          opt.ulim, opt.sp, ...
                         [opt.wb opt.wh],opt.N,opt.type,opt.cancelzero);
                 
       [y, u, t, r] = fpid_optimize_sim([Kp,Ki,lam,Kd,mu], ...
                                   Z_p, sopt);
       
       % Compute error
       e  = r - y;
       
       % Performance metric calculation
       switch (lower(opt.optim))
            case 'ise'
                z=trapz(t, e.^2);
            case 'iae'
                z=trapz(t,abs(e));
            case 'itse'
                z=trapz(t,(t.*(e.^2)));
            case 'itae'
                z=trapz(t,t.*abs(e));
            otherwise
                error('Unknown performance metric!');
       end
       
       % If control law metric is to be used, compute it here
       if ~isempty(opt.wgt)
            
           u_avg   = 1/length(u) * sum(u);
           u_sumsq = sum((u-u_avg).^2);
           
           % Compute the final weighted objective function
           z       = (1-opt.wgt) * indexWeight * z + ...
                      opt.wgt / length(u) * u_sumsq;
       end
       
    else
        % Otherwise use linear simulation
        if fleq(opt.sp, 1)
            e = 1 - step(feedback(Z_p*Z_c,1), t);
        else
            u = opt.sp * ones(length(t), 1);
            e = opt.sp - lsim(feedback(Z_p*Z_c,1), u, t);
        end
        
        % Performance metric calculation
        switch (lower(opt.optim))
            case 'ise'
                z=e'*e*dt;
            case 'iae'
                z=sum(abs(e)*dt);
            case 'itse'
                z=(t.*e'*dt)*e;
            case 'itae'
                z=sum(t'.*abs(e)*dt);
            otherwise
                error('Unknown performance metric!');
        end
        
        % Performance index weighting factor
        z = indexWeight * z;
        
    end
   
end

