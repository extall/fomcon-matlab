function [R,C,L,K,results] = frac_imp_rc_foster2_abgen(params)
%FRAC_IMP_RC_ABGEN_FOSTER2 Generate component values for Foster II RC imp.

    str_fun  = params.structure;
    orig_fun = params.model;
    R1       = params.R1;
    C1       = params.C1;
    varphi   = params.varphi;
    N        = params.N;
    
    % Fetch FOTF parameters
    [fa, fna, fb, fnb] = fotfparam(orig_fun);
    
    if (numel(fa)  ~= 1  || ...
        numel(fna) ~= 1  || ...
        numel(fb)  ~= 1  || ...
        numel(fnb) ~= 1)
              
       % Unsupported type of system
       error('FracImpRcAbGenFoster2:UnsupportedSystem', ...
             'Unsupported system given.');
      
    end

    % Operator order
    if fleq(fnb,0)
        alpha = fna;
    elseif fleq(fna,0)
        alpha = -fnb;
    else
        % Unsupported type of system
       error('FracImpRcAbGenFoster2:UnsupportedSystem', ...
             'Unsupported system given.');
    end

    % Ratio (related to maximum phase ripple)
    r_ab = 0.24 / (1+varphi);
    
    % Generation parameters
    a = 10^(alpha*log10(r_ab));
    b = r_ab / a;
    
    % Component value generation
    R = R1 .* a.^(0:1:N-1)';
    C = C1 .* b.^(0:1:N-1)';
    
    % Corection components
    Rp = R1*(1-a)/a; Cp = C1*b^N/(1-b);
    
    % Return component values
    R = [Rp; R]; C = [Cp; C];
    
    % Get structure function
    struct_fun = str2func(str_fun);
    
    % Compute gain
    Z    = struct_fun(R,C,[],1);
    w_av = 1/(R1*C1*(r_ab)^(floor(N/2)-1))*sqrt(a);
    
    K = abs(squeeze(freqresp(orig_fun, w_av))) / ...
        abs(squeeze(freqresp(Z, w_av)));
    
    % Unset parameters
    L       = [];
    results = [];
    
end

