function [R,C,K,results] = frac_imp_rc_optimize_abRC_foster2(params)
%FRAC_IMP_RC_OPTIMIZE Obtain R and C component values for an analog
%  approximation of a given fractional-order system

    % Get optimization parameters
    str_fun  = params.structure;     % Used ciruit structure
    orig_fun = params.G;             % Original fractional-order system
    w        = params.w;             % Fequency response points
    N        = params.N;             % Number of components
    
    % Results array
    results = {};
    
    % Optimization parameters
    op = set_options('Display', 'on');
    
    % Begin optimization
    x = GODLIKE(@(x) obj_func(x, str_fun, orig_fun, w, N), ...
                12, ...
                [0 0 1   1e-12 1e-6], ...
                [1 1 1e6 1e-6  1e+6], ...
                {'GA'}, ...
                op);
          
    % Set result
    [R0, C0, R, C, K] = get_rc_values(x,N);
    R = [R0; R];
    C = [C0; C];

end

function [R0, C0, R, C, K] = get_rc_values(x, N)

    % Generation parameters
    a = x(1); b = x(2);
    R1 = x(3); C1 = x(4);
    K = x(5);
    
    % Component values
    R = R1*a.^(0:N-1)';
    C = C1*b.^(0:N-1)';
    R0 = R1*(1-a)/a;
    C0 = C1*b^N/(1-b);

end

function f = obj_func(x, str_fun, orig_fun, w, N)

    % Get structure function
    struct_fun = str2func(str_fun);
    
    % Build the RC network
    [R0, C0, R, C, K] = get_rc_values(x,N);   % Get R,C network values
    Z = struct_fun([R0; R],[C0; C])*K;        % Get impedance K Z(s)
    
    % Compute the objective function value from the frequency response
    r_orig = squeeze(freqresp(orig_fun, w))';
    r_this = squeeze(freqresp(Z, w));
    
    % Index calculation
    f = sum(abs(r_orig - r_this).^2)/length(w);
    
    
end

