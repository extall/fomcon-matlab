function z = fracidfun(x, y, u, t, opt)
%FRACIDFUN Compute fractional time-domain identification cost function.
%
%   Usage: z = fracidfun(x,y,u,t,opt)
%
%          opt.type - 'gl', 'oust' or 'ref'
%          opt.finddelay - true/false, whether to use the delay parameter
%                          for model estimation
%          opt.wb, opt.wh, opt.N - parameters to be used
%                                  with 'oust' or 'ref'
%          opt.fix - identification type, 'n' for free identification,
%                    'c' to fix coefficients, 'e' to fix exponents
%          opt.fixpoly - a matrix [zeroPolyLength; polePolyLength],
%                        if length is zero, the polynomial is considered
%                        fixed and derived from initial model. Both lengths
%                        cannot be zero at once.
%          opt.G - initial model
%          opt.alg - identification algorithm (1: TRR, 2: LM)

    % Check fixpoly
    bSize = opt.fixpoly(1);
    aSize = opt.fixpoly(2);
    
    % Cannot identify - both polynomials fixed
    if bSize == 0 && aSize == 0
        error('FRACIDFUN:BothPolynomialsFixed', ...
              'Cannot identify model because both polynomials are set to be fixed');
    end
    
    % Get initial model parameters
    [ia, ina, ib, inb] = fotfparam(opt.G);

    % NB! If LM algorithm is used, the variables are transformed,
    % so must be transformed back by means of x = z^2, use the
    % corresponding function for that at the end of this file.
    
    % Extract delay from x if required
    if opt.finddelay
        ioDelay = x(1);
        ioDelay = convertFromLM(ioDelay, opt.alg);
        x = x(2:end);
    else
        ioDelay = 0;
    end
    
    % Extract static gain
    if opt.findstaticgain
        % Extract the gain from the optimization vector
        dcGain = x(1);
        x = x(2:end);
        
        % Remove the original free terms
        ia  = ia(1:end-1);
        ina = ina(1:end-1);
        ib  = ib(1:end-1);
        inb = inb(1:end-1);
    else
        dcGain = 1;
    end

    % Pole polynomial is fixed
    if bSize > 0 && aSize == 0
        
        [b, nb] = fracidfun_getpolyparam(opt.fix, bSize, x, ib, inb);
        % Only LM-convert unless exponents are fixed
        if ~strcmpi(opt.fix, 'e'), nb = convertFromLM(nb, opt.alg); end
        a = ia;
        na = ina;
        
    % Zero polynomial is fixed
    elseif bSize == 0 && aSize > 0
        
        b = ib;
        nb = inb;
        [a, na] = fracidfun_getpolyparam(opt.fix, aSize, x, ia, ina);
        if ~strcmpi(opt.fix, 'e'), na = convertFromLM(na, opt.alg); end
        
    % Free identification
    elseif bSize > 0 && aSize > 0
        
        [a, na, b, nb] = fracidfun_getfotfparam(opt.fix, bSize, aSize, x, ia, ina, ib, inb);
        if ~strcmpi(opt.fix, 'e')
            na = convertFromLM(na, opt.alg);
            nb = convertFromLM(nb, opt.alg);
        end
        
    end
    
    % Add the static gain, if required
    if opt.findstaticgain
        b = dcGain*[b 1];
        a = [a 1];
        nb = [nb 0];
        na = [na 0];
    end
    
    % Get identified fotf object
    G = fotf(a,na,b,nb,ioDelay);

    % Build model based on type
    switch lower(opt.type)
        case 'gl'
            y_id = lsim(G,u,t);
        case 'oust'
            G = oustapp(G,opt.wb,opt.wh,opt.N,'oust');
            
            % Convert to proper model for simulation
            G = toproper(G, opt.wh);
            
            y_id = lsim(G,u,t);
        case 'ref'
            G = oustapp(G,opt.wb,opt.wh,opt.N,'ref');
            
            % Convert to proper model for simulation
            G = toproper(G, opt.wh);
            
            y_id = lsim(G,u,t);
        otherwise
            error('Unknown simulation type specified!');
    end
    
    % Return error
    z = y - y_id;
    
end

% Returns polynomial parameters based on desired identification type
function [p, np] = fracidfun_getpolyparam(fix, vec_len, vec, ip, inp)

    switch(lower(fix))
        case 'n'
            % Free identification
            p = vec(1:vec_len/2);
            np = vec(vec_len/2+1:end);
        case 'e'
            % Fix exponents
            p = vec;
            np = inp;
        case 'c'
            % Fix coefficients
            p = ip;
            np = vec;
    end

end

% Returns both polynomial parameters based on desired identification type
function [a, na, b, nb] = fracidfun_getfotfparam(fix, bSize, aSize, vec, ia, ina, ib, inb)
    
    switch(lower(fix))
        case 'n'
            b = vec(1:bSize/2);
            nb = vec(bSize/2+1:bSize);
            a = vec(bSize+1:aSize/2+bSize);
            na = vec(aSize/2+bSize+1:end);
        case 'e'
            b = vec(1:bSize);
            nb = inb;
            a = vec(bSize+1:end);
            na = ina;
        case 'c'
            b = ib;
            nb = vec(1:bSize);
            a = ia;
            na = vec(bSize+1:end);
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