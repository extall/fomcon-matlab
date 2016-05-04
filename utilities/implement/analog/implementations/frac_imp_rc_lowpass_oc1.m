function [R,C,L,K,results] = frac_imp_rc_lowpass_oc1(params)
%FRAC_IMP_RC_LOWPASS_OC1 Generate component values for Cauer I RC LP imp.

    str_fun  = params.structure;
    orig_fun = params.model;
    R1       = params.R1;
    N        = params.N;
    w        = params.w;
    
    % Fetch FOTF parameters
    [fa, fna, fb, fnb] = fotfparam(orig_fun);
    
    if (numel(fa)  ~= 2  || ...
        numel(fna) ~= 2  || ...
        numel(fb)  ~= 1  || ...
        numel(fnb) ~= 1  || ...
        ~fleq(fnb,0)        ...
       )
              
       % Unsupported type of system
       error('FracImpRcAbGenFoster2:UnsupportedSystem', ...
             'Unsupported system given.');
      
    end
    
    % Fetch fractional operator order and RC constant
    alpha        = fna(1);
    scale_factor = fa(1) / R1;
    
    % Obtain Oustaloup approximation
    G = 1/zpk((scale_factor*oustafod(alpha,N,w(1),w(end))));
    
    % Check orders of zero/pole polynomials
    if length(zpk(G).p{1}) ~= length(zpk(G).z{1})
        error('FracImpRcLowpassOc1:BadApproximation', ...
              'Approximation failed.');
    end
    
    % Develop approximation into CFE
    Rt = []; Ct = []; R = {}; C = {};
    q = polycfe(G);
    for n=1:length(q)
        
        tt = q{n};
        if length(tt) > 1       % Capacitor
            Ct(end+1) = tt(1);
        elseif length(tt) == 1  % Resistor
            Rt(end+1) = tt;
        end
            
    end
    
    % Assign resistor and capacitor arrays
    R{1} = R1;
    R{2} = Rt;
    C{1} = [];
    C{2} = Ct;
    
    % Compute gain
    struct_fun = str2func(str_fun);
    Z          = struct_fun(R,C,[],{1,1});
    w_av       = sqrt(w(1)*w(end));
    Kg         = abs(squeeze(freqresp(orig_fun, w_av))) / ...
                 abs(squeeze(freqresp(Z, w_av)));
    
    K{1} = Kg;
    K{2} = 1;
             
    L = [];
    results = [];
    
end

