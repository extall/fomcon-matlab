function [R,C,results] = frac_imp_rc_optimize_ab_foster2(params)
%FRAC_IMP_RC_OPTIMIZE Obtain R and C component values for an analog
%  approximation of a given fractional-order system

    % Get optimization parameters
    str_fun  = params.structure;     % Used ciruit structure
    orig_fun = params.G;             % Original fractional-order system
    w        = params.w;             % Fequency response points
    a        = params.a;             % a, b - geometric sequence parameters
    b        = params.b;
    R1       = params.R1;            % Base resistor and capacitor values
    C1       = params.C1;
    N        = params.N;             % Number of components (total number: 2N+2)
    
    % Results array
    results = {};
    
    % Initial values
    init  = [a b R1 C1];
    
    % Optimization parameters
    op = optimset('Display', 'iter');
    
    % Begin optimization
    x = optimize(@(x) obj_func(x, str_fun, orig_fun, w, N), ...
                 init, ...
                 [0 0 1 1e-12], ...
                 [1 1 1e6 1e-9], ...
                 [], [], [], [], [], [], ...
                 op ...
              );
          
    % Set result
    [R0, C0, R, C] = get_rc_values(x(1),x(2),x(3),x(4),N);
    R = [R0; R];
    C = [C0; C];
    
    a = x(1)
    b = x(2)
    R1 = x(3)
    C1 = x(4)

end

function [R0, C0, R, C] = get_rc_values(a, b, R1, C1, N)

    R = R1*a.^(0:N-1)';
    C = C1*b.^(0:N-1)';
    R0 = R1*(1-a)/a;
    C0 = C1*b^N/(1-b);

end

function f = obj_func(x, str_fun, orig_fun, w, N)

    % Get structure function
    struct_fun = str2func(str_fun);
    
    % Build the RC network
    a = x(1); b=x(2); R1=x(3); C1=x(4);
    [R0, C0, R, C] = get_rc_values(a,b,R1,C1,N);   % Get R,C network values
    Z = struct_fun([R0; R],[C0; C]);               % Get impedance Z(s)
    
    % Compute the objective function value from the frequency response
    r_orig = squeeze(freqresp(orig_fun, w))';
    r_this = squeeze(freqresp(Z, w));
    
    % Index calculation
    f = sum(abs(r_orig - r_this).^2)/numel(w);
    
end

function [c, ceq] = nonlin_con(x, orig_fun, w, R1, C1, N)

end

