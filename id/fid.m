function [a,na,b,nb,L,Gid] = fid(fs,gparam,idd,npoints,type,fixpoly,limits,op)
%FID Fractional-order Transfer Function SISO system identification
%   
%Usage: [A,NA,B,NB,L,GID] = FID(FSIM|G,GPARAM,IDD,NPOINTS,TYPE,FIXPOLY,LIMITS,OP)
%
%where  A,NA,B,NB, L and GID - identified system coefficients, input/output
%                              delay from delay term exp(-Ls) and FOTF object.
%                              If the algorithm is 
%
%       FSIM | G    - either FSPARAM structure to use approximations for
%                     identification with the initial function or the
%                     initial function as a FOTF object. If the FOTF object
%                     has a set ioDelay value, it will also be identified
%                     in the ensuing optimization. To ignore the delay
%                     parameter set it to [] (i.e. empty matrix) in the initial
%                     model.
%       GPARAM      - a cell array with explicit model parameters: {K, L},
%                     where K is the static gain and L is the input-output
%                     delay in seconds. Specifying this input argument as
%                     [] (or each individual parameter, e.g. {K, []}) will
%                     remove the parameter(s) from the identification problem.
%                     NOTE! The static gain will ONLY be identified in case
%                     of free identification (all parameters of both
%                     polynomials are identified).
%       IDD         - a FIDATA sturcture with collected system samples.
%       ----------------------------------------------------------------
%       Optional arguments
%       ----------------------------------------------------------------
%       NPOINTS     - number of points to use for identification; this
%                     allows to reduce the time it takes to estimate a
%                     model, i.e. for obtaining an initial estimate, at an
%                     expense of accuracy. If the number of points
%                     specified here is larger than half the number of
%                     points in the dataset, all available points will be
%                     used for identification. Specifying 0 or [] will
%                     yield the same effect, that is all points are used,
%                     and it is also the default behavior.
%       TYPE        - identification type: 'n' - free identification,
%                     'c' - fix coefficients, 'e' - fix exponents
%                     (default: 'n')
%       FIXPOLY     - a vector with two values: [BFIX; AFIX], where BFIX
%                     and AFIX can be 1, in which case the corresponding
%                     polynomial is fixed during identification, or 0.
%                     Note that in case BFIX = AFIX = 1 the initial model
%                     will be immediately returned with no identification
%                     conducted. Default: [0; 0]
%       LIMITS      - a cell array with two vectors, which may also be
%                     empty, containing polynomial coefficient and 
%                     exponent limits in the form {[CMIN;CMAX],[EMIN;EMAX]}.
%                     Default: [] (i.e. no limits are imposed.)
%       OP          - additional optimization options for lsqnonlin
%                     (use optimset). Note, that you can set the preferred
%                     optimization algorithm here. This is done as follows:
%                     OP.IdentificationAlgorithm = 'trr':
%                     Trust-Region-Reflective algorithm is used;
%                     OP.IdentificationAlgorithm = 'lm':
%                     Levenberg-Marquardt algorithm is used, for the latter
%                     OP.Lambda determines the lambda parameter, 0.01 by
%                     default. NB! The 'lm' algorithm does not handle bound
%                     constraints! The upper/lower limits of search
%                     parameter space will be discarded!
%
%       See also: fsparam, fidata, optimset, lsqnonlin

    % Load configuration parameters
    config = fomcon('config');
    numSigDig = config.Core.General.Model_significant_digits;
	epsi   = config.Core.General.Internal_computation_accuracy;

    % Check simulation parameters
    if isa(fs, 'fsparam')
        G = fs.plant;
        opt.type = fs.approx;
        opt.wb = fs.w(1);
        opt.wh = fs.w(2);
        opt.N = fs.N;
    elseif isa(fs, 'fotf')
        G = fs;
        opt.type = 'gl';
    else
        error('FID:WrongSimulationParameters', ...
              'Unexpected simulation object type!');
    end
    
    % Check identification array
    if ~isa(idd, 'fidata') || ...
        size(idd.y,1) ~= size(idd.u,1) || ...
        size(idd.u,1) ~= size(idd.t,1)
            error('FID:WrongIdentificationData', ...
                'Bad identification data or sizing error!');
    else
        y  = idd.y;
        u  = idd.u;
        t  = idd.t;
        dt = idd.dt;
    end
    
    % Check value of the response
    if abs(y(1)) > eps
       warning('FID:NONZERORESPONSE', ...
           'Nonzero first element in response dataset.');
    end
    
    % Check optional arguments
    
    % Optimization options and algorithm choice
    % 1: Trust-Region-Reflective
    % 2: Levenberg-Marquardt
    alg = 1;
    if nargin < 8
        op = optimset('Display','iter');
    else
        % Determine the algorithm and set parameters accordingly
        if cfieldexists(op, 'IdentificationAlgorithm')
           switch lower(op.IdentificationAlgorithm)
               case 'lm'
                   alg = 2;
                   op.Algorithm = 'levenberg-marquardt';
                   if cfieldexists(op, 'Lambda') && isnumeric(op.Lambda)
                      op.Algorithm = {'levenberg-marquardt', op.Lambda}; 
                   end
           end
        end
    end
    
    opt.alg = alg;
    
    % Coefficient/exponent limits
    if nargin < 7
       limits = []; 
    end
    
    % Polynomial fix
    if nargin < 6 || isempty(fixpoly)
        fixpoly = [0; 0];
    end
    
    % Identification type
    if nargin < 5 || isempty(type)
        type = 'n';
    end
    
    % Identification type
    if nargin < 4
        npoints = [];
    end
    
    % Everything else supplied?
    if nargin < 3
        error('FID:NotEnoughInputArguments', ...
              'Not enough input arguments.');
    end
    
    % Separate model parameters K and L, and limits
    if ~isempty(gparam)
        optK = gparam{1};
        optL = gparam{2};
    end
    
    if ~isempty(limits)
        clim = limits{1};
        elim = limits{2};
    else
        clim = [];
        elim = [];
    end
    
    % Check polynomial fix options
    if fixpoly(1) == 1 && fixpoly(2) == 1
       % Both fixed, nothing to identify
       [a, na, b, nb] = fotfparam(G);
       Gid = G;
       return;
    end
    
    % Check exponent limits
    if ~isempty(elim)
        if elim(1) > elim(2)
            
            temp = elim(1);
            elim(1) = elim(2);
            elim(2) = temp;
            
            warning('FID:ExponentLimitsSwapped', ...
                    'Swapping exponent min and max limits');
        end 
    end
    
    % Check coefficient limits
    if ~isempty(clim)
        if clim(1) > clim(2)
           
            temp = clim(1);
            clim(1) = clim(2);
            clim(2) = temp;
            
            warning('FID:CoefficientLimitsSwapped', ...
                    'Swapping coefficient min and max limits');
        end
    end
    
    % Check the number of used points
    if ~isempty(npoints) && npoints ~= 0 && npoints < numel(t)/2
        % Create new array with specified inter-point distance
        interpoints = floor(numel(t)/npoints);
        t = t(1:interpoints:end);
        y = y(1:interpoints:end);
        u = u(1:interpoints:end);
    end
    
    % Initial model
    opt.G = G;
    
    % Fix type
    opt.fix = type;
    
    % Delay parameter estimation
    if ~isempty(optL)
       opt.finddelay = 1;
       ioDelay_l = 0;
       ioDelay_u = t(end)-dt;       % Maximum delay limited by time vector
    else
       opt.finddelay = 0;
       ioDelay_l = [];
       ioDelay_u = [];
    end
    
    % DC gain estimation
    if ~isempty(optK)
        opt.findstaticgain = 1;
        dcGain_l = epsi/2;
        dcGain_h = Inf;
    else
        opt.findstaticgain = 0;
        dcGain_l = [];
        dcGain_h = [];
    end
    
    % If DC gain is to be estimated, only free
    % parameter identification is allowed
    if opt.findstaticgain && ...
       ((fixpoly(1) == 1 || fixpoly(2) == 1) || ...
       ~strcmpi(type,'n'))
   
            error('FID:StaticGainEstimationOnFixedModel', ...
                  'FID() cannot estimate static gain K for fixed polynomial parameters');
    end
       
    % Determine exponent bounds
    if isempty(elim)
       leb = 0;
       ueb = Inf;
    else
       leb = elim(1);
       ueb = elim(2);
    end
    
    % Check lower exponent bound
    if leb < 0
        leb = 0;
        warning('FID:BadExponentLowerLimit', ...
                'Invalid lower limit set for exponents; setting to zero');
    end
    
    % Determine coefficient bounds
    if isempty(clim)
       lcb = -Inf;
       ucb = Inf;
    else
       lcb = clim(1);
       ucb = clim(2);
    end
    
    % Initial guess parameters
    [xa, xna, xb, xnb] = fotfparam(G);
    
    % Remove free term places in case of
    % a specified static gain parameter
    if opt.findstaticgain
        xa  = xa(1:end-1);
        xna = xna(1:end-1);
        xb  = xb(1:end-1);
        xnb = xnb(1:end-1);
    end  
    
    % If LM algorithm is used, we must use a variable transformation
    % method---providing this way a lower bound---so that exponents are
    % always positive. The transformation is x = z^2, so for initial
    % estimates we do z = sqrt(x).
    if (alg==2)
       xna = sqrt(xna);
       xnb = sqrt(xnb);
    end
    
    % Set initial guess, lower and upper bounds
    if (fixpoly(1) == 0 && fixpoly(2) == 1)
        % Pole polynomial is fixed
        
        % Element fix choices
        switch(lower(type))
            case 'n'
                % Free identification
                
                % Initial guess
                x0 = [xb xnb];
                
                % Limits
                lb = [lcb*ones(1,length(x0)/2) leb*ones(1,length(x0)/2)];
                ub = [ucb*ones(1,length(x0)/2) ueb*ones(1,length(x0)/2)];
            case 'e'
                % Fix exponents
                x0 = [xb];
                lb = [lcb*ones(1,length(x0))];
                ub = [ucb*ones(1,length(x0))];
            case 'c'
                % Fix coefficients
                x0 = [xnb];
                lb = [leb*ones(1,length(x0))];
                ub = [ueb*ones(1,length(x0))];
        end
        
        % Set sizes
        opt.fixpoly = [length(x0);0];
        
    elseif (fixpoly(1) == 1 && fixpoly(2) == 0)
        % Zero polynomial is fixed
        switch(lower(type))
            case 'n'
                x0 = [xa xna];
                lb = [lcb*ones(1,length(x0)/2) leb*ones(1,length(x0)/2)];
                ub = [ucb*ones(1,length(x0)/2) ueb*ones(1,length(x0)/2)];
            case 'e'
                x0 = [xa];
                lb = [lcb*ones(1,length(x0))];
                ub = [ucb*ones(1,length(x0))];
            case 'c'
                x0 = [xna];
                lb = [leb*ones(1,length(x0))];
                ub = [ueb*ones(1,length(x0))];
        end
        
        opt.fixpoly = [0;length(x0)];
        
    else
        % Full model identification
        switch(lower(type))
            case 'n'
                x01 = [xb xnb];
                x02 = [xa xna];
                x0 = [x01 x02];
                lb = [lcb*ones(1,length(x01)/2) leb*ones(1,length(x01)/2) ...
                      lcb*ones(1,length(x02)/2) leb*ones(1,length(x02)/2)];
                ub = [ucb*ones(1,length(x01)/2) ueb*ones(1,length(x01)/2) ...
                      ucb*ones(1,length(x02)/2) ueb*ones(1,length(x02)/2)];
            case 'e'
                x01 = xb;
                x02 = xa;
                x0 = [x01 x02];
                lb = [lcb*ones(1,length(x0))];
                ub = [ucb*ones(1,length(x0))];
            case 'c'
                x01 = xnb;
                x02 = xna;
                x0 = [x01 x02];
                lb = [leb*ones(1,length(x0)-1) 1e-10];
                ub = [ueb*ones(1,length(x0)-1) 1e-9];
        end
        
        opt.fixpoly = [length(x01); length(x02)];
         
    end
    
    % Add static gain, if present
    if ~isempty(optK)
        x0 = [optK x0];
        lb = [dcGain_l lb];
        ub = [dcGain_h ub];
    end
    
    % Add delay, if present
    if ~isempty(optL)
        % In case of the LM algorithm, the delay must be nonnegative
        optL = convertToLM(optL, alg);
        
        x0 = [optL x0];
        lb = [ioDelay_l lb];
        ub = [ioDelay_u ub];
    end
    
    % Disable pbar if enabled
    sp = show_pbar(); if sp, show_pbar('off'); end
    
    % Measure time
    elTime = tic;
    
    % Run the identification algorithm
    switch alg
        case 1
            [x_opt, resNorm, res] = lsqnonlin(@(x) fracidfun(x,y,u,t,opt), ...
                                                x0, ...
                                                lb, ...
                                                ub, ...
                                                op); 
        case 2
            [x_opt, resNorm, res] = lsqnonlin(@(x) fracidfun(x,y,u,t,opt), ...
                                                x0, ...
                                                [], ...
                                                [], ...
                                                op);
    end

    % Make sure that we have a horizontal vector
    if size(x_opt,1)>size(x_opt,2) x_opt = x_opt'; end
    
    % Display elapsed time
    elTime = stohms(toc(elTime));
    elTime = ['Elapsed time: ' elTime];
    disp(elTime);
    disp(char('-'*ones(1,length(elTime))));
    
    disp(['Residual norm: ' num2str(resNorm,numSigDig)]);
    disp('Identification completed.');
    
    % Re-enable pbar
    if sp, show_pbar('on'); end
    
    % Set initial model parameters by default
    [a, na, b, nb, L] = fotfparam(G);
    
    % Fetch delay
    if opt.finddelay
        L = x_opt(1);
        x_opt = x_opt(2:end);
        
        % Convert delay in case of LM algorithm
        L = convertFromLM(L, alg);
        
	else
		L = 0;
    end
    
    % Fetch static gain
    if opt.findstaticgain
        K = x_opt(1);
        x_opt = x_opt(2:end);
    else
        K = 1;
    end
    
    % And update those which were optimized
    if (fixpoly(1) == 0 && fixpoly(2) == 1)
        % Pole polynomial is fixed
        switch(lower(type))
            case 'n'
                b = x_opt(1:length(x_opt)/2);
                nb = x_opt(length(x_opt)/2+1:end);
                nb = convertFromLM(nb,alg);
            case 'e'
                b = x_opt;
            case 'c'
                nb = x_opt;
                nb = convertFromLM(nb,alg);
        end
        
    elseif (fixpoly(1) == 1 && fixpoly(2) == 0)
        % Zero polynomial is fixed
        switch(lower(type))
            case 'n'
                a = x_opt(1:length(x_opt)/2);
                na = x_opt(length(x_opt)/2+1:end);
                na = convertFromLM(na,alg);
            case 'e'
                a = x_opt;
            case 'c'
                na = x_opt;
                na = convertFromLM(na,alg);
        end
        
    else
        % Full model identification
        switch(lower(type))
            case 'n'
                b = x_opt(1:opt.fixpoly(1)/2);
                nb = x_opt(opt.fixpoly(1)/2+1:opt.fixpoly(1));
                a = x_opt(opt.fixpoly(1)+1:opt.fixpoly(1)+opt.fixpoly(2)/2);
                na = x_opt(opt.fixpoly(1)+opt.fixpoly(2)/2+1:end);
                na = convertFromLM(na,alg);
                nb = convertFromLM(nb,alg);
            case 'e'
                b = x_opt(1:opt.fixpoly(1));
                a = x_opt(opt.fixpoly(1)+1:end);
            case 'c'
                nb = x_opt(1:opt.fixpoly(1));
                na = x_opt(opt.fixpoly(1)+1:end);
                na = convertFromLM(na,alg);
                nb = convertFromLM(nb,alg);
        end
    end
    
    % Return the model with static gain
    if opt.findstaticgain
       b = K*[b 1];
       a = [a 1];
       nb = [nb 0];
       na = [na 0];
    end
    
    % Return fotf
    Gid = fotf(a,na,b,nb,L);
    
end

% Conversion TO LM coordinates (lower bounds)
function nq1 = convertToLM(nq, alg)
    if (alg==2)
        nq1=sqrt(nq);
    else
        nq1=nq;
    end 
end

% Conversion FROM LM coordinates (lower bounds)
function nq1 = convertFromLM(nq, alg)
    if (alg==2)
        nq1=nq.^2;
    else
        nq1=nq;
    end
end

