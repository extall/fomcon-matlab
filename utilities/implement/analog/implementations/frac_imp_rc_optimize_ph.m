function [R,C,K,results] = frac_imp_rc_optimize_ph(params)
%FRAC_IMP_RC_OPTIMIZE Obtain R and C component values for an analog
%  approximation of a given fractional-order system

    % Get optimization parameters
    str_fun  = params.structure;     % Used ciruit structure
    orig_fun = params.G;             % Original fractional-order system
    w        = params.w;             % Fequency response points
    N        = params.N;             % Number of discrete components
    varphi   = params.varphi;        % Allowed phase ripple [deg]
    
    % Results array
    results = {};
    
    % Optimization parameters
    op = set_options('display', 'on');
    
    % TEST
    orig_fun=oustapp(orig_fun,w(1),w(end),10,'oust');
    r_orig = squeeze(freqresp(orig_fun, w));
    % END TEST
    
    % Begin optimization
    x = GODLIKE(@(x) obj_func(x, str_fun, r_orig, w, N, varphi), ...
                 12, ...
                 [1*ones(1,N)    1e-12*ones(1,N)], ...
                 [1e6*ones(1,N) 10e-6*ones(1,N)], ...
                 {'GA'}, ...
                 op ...
              );
          
    % Set result
    [R, C] = get_rc_values(x,N);
    
    % Get structure function
    struct_fun = str2func(str_fun);
    
    Z = struct_fun(R,C);
    w_av = sqrt(w(1)*w(end));
    K = abs(squeeze(freqresp(orig_fun, w_av))) / ...
        abs(squeeze(freqresp(Z, w_av)));
    
end

function [R,C] = get_rc_values(x,N)

    R  = x(1:N);
    C  = x(N+1:end);
    
end

function f = obj_func(x, str_fun, r_orig, w, N, varphi)

    % Get structure function
    struct_fun = str2func(str_fun);
    
    % Build the RC network
    [R,C] = get_rc_values(x,N);      % Get R,C network values
    
    % Find correcting coefficient
    Z = struct_fun(R,C);
%     w_av = sqrt(w(1)*w(end));
%     K = abs(squeeze(freqresp(orig_fun, w_av))) / ...
%         abs(squeeze(freqresp(Z, w_av)));
%     
%     Z = Z * K;

    % Try different approach: constant phase
    % r_orig = squeeze(freqresp(orig_fun, w));
    r_this = squeeze(freqresp(Z, w));
     
    % Cost function as square norm of phase angle difference within
    % a tolerance set by the varphi parameter
    e_abs = rad2deg(abs(angle(r_orig) - angle(r_this)));
    %e_phi(e_abs>varphi)=e_abs(e_abs>varphi);
    % f = sum(e_phi);
    f = sum(e_abs.^2);

end

