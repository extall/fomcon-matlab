function [R,C,K,results] = frac_imp_rc_optimize_foster2(params)
%FRAC_IMP_RC_OPTIMIZE Obtain R and C component values for an analog
%  approximation of a given fractional-order system

    % Get optimization parameters
    str_fun  = params.structure;     % Used ciruit structure
    orig_fun = params.G;             % Original fractional-order system
    w        = params.w;             % Fequency response points
    N        = params.N;             % Number of discrete components
    
    % Results array
    results = {};
    
    % Band-limit the original function (i.e. create approximation)
    orig_fun = oustapp(orig_fun,w(1),w(end),N,'oust');
    
    % Optimization parameters
    op = set_options('display', 'on', 'Coding', 'real');
    
    % Begin optimization
    x = GODLIKE(@(x) obj_func(x, str_fun, orig_fun, w, N), ...
                 12, ...
                 [1e-6 0   0      1*ones(1,N)   1e-12*ones(1,N)], ...
                 [+1e6 1e6 10e-6  1e6*ones(1,N) 10e-6*ones(1,N)], ...
                 {'GA'}, ...
                 op ...
              );
          
    % Set result
    [K, Rp, Cp, R, C] = get_rc_values(x,N);
    R = [Rp R]; C = [Cp C];
    
end

function [K,Rp,Cp,R,C] = get_rc_values(x,N)
    
    K  = x(1); Rp = x(2); Cp = x(3); x = x(4:end);
    R  = x(1:N);
    C  = x(N+1:end);
    
end

function f = obj_func(x, str_fun, orig_fun, w, N)

    % Get structure function
    struct_fun = str2func(str_fun);
    
    % Build the RC network
    [K,Rp,Cp,R,C] = get_rc_values(x,N);          % Get R,C network values
    Z = struct_fun([Rp R],[Cp C])*K;
    
%     % Compute the cost function value from the frequency response
%     r_orig = squeeze(freqresp(orig_fun, w))';
%     r_this = squeeze(freqresp(Z, w));
%     
%     % Index calculation
%     f = sum(abs(r_orig - r_this).^2)/length(w);

      % Try different approach: constant phase
      r_orig = squeeze(freqresp(orig_fun, w));
      r_this = squeeze(freqresp(Z, w));
      
      % Sum of magnitude difference squares
      f1 = sum((abs(r_orig)-abs(r_this)).^2);
      
      % Sum of phase difference squares
      f2 = sum((angle(r_orig) - angle(r_this)).^2);
      
      % Full cost function is a sum (weighted?)
      f = f1 + f2;

end

