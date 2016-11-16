function [G, J, err] = fpzp_new(b, a, alpha, N, w, retol)
%FPZP_NEW Fractional power zero-pole pair approximation by Newton method
%
%   Approximates a fractional power pole-zero pair in form 
%
%                                      alpha
%                            (b*s + 1)
%                            ---------
%                            (a*s + 1)
%
%  where alpha is the fractional order (power), b and a are coefficients
%  such that w_zero = 1/b, w_pole = 1/a.
%
%  Usage: [G, J, err] = fpzp_new(b, a, alpha, N, w, retol)
%
%         where optional parameters are
%
%                   N (default: 2) - approximation order,
%
%                   w - [wb; wh] - frequency range, where index J is
%                                  calculated. Defaults to [0.0001; 10000].
%                   
%                   retol (default: []) - tolerance for model reduction
%                                         using a matching zero-pole
%                                         reduction method realized by
%                                         applying the minreal() function.
%
%  See also: carlapp, minreal

    % The following sets the maximum power M for ALPHA = 1/M,
    % default is 10, and while it provides limited resolution, taking
    % larger values without proper model simplification may result in a
    % system of a very high order.
    M = 10;

    % Check input argument number
    if nargin < 3
        error('FPZPNEW:NotEnoughInputArguments', ...
              'Not enough input arguments.');
    end

    % Check input arguments
    if nargin < 4, N = 2; end                         % Approximation order
    
    if nargin < 5 || isempty(w) || max(size(w)) ~= 2  % Frequency range
        w = [0.0001; 10000];
    end
    
    if nargin < 6 || isempty(retol), retol = []; end  % Reduction tolerance
    
    % Check coefficients
    if max(size(b)) ~= 1
        warning('FPZPNEW:BMUSTBESCALAR', ...
                'Coefficient B must be a scalar, using first cell value.');     
        b = b(1);
    end
    
    if max(size(a)) ~= 1
        warning('FPZPNEW:AMUSTBESCALAR', ...
                'Coefficient A must be a scalar, using first cell value.');     
        a = a(1);
    end
    
    % Model
    Gi = (b * zpk('s') + 1) / ...
         (a * zpk('s') + 1);
    
    % Save initial value of alpha and Gi
    alpha_init = alpha;
    Gi_init = Gi;
     
    % Check alpha
    if alpha < 0
        Gi = inv(Gi);
        alpha = -alpha;
    end
    
    % Extract integer order
    if alpha >= 1
        v = fix_s(alpha);
        alpha = alpha - v;
    else
        v = 0;
    end
    
    % If only integer order is present, return result immediately
    if fleq(alpha, 0) && ~fleq(v, 0)
        G   = Gi^v;
        J   = 0;
        err = 0;
        return;
    end
    
    % Initialize system
    G = 1;
    
    % Check alpha resolution
    if fleq(alpha, fix_s(alpha*10)/10)
        % Use optimized algorithm for accuracy: use expansion by 1/2,
        % 1/5 and 1/10
        
        % 0.2 multiplier
        v02 = 0;
        
        % Approximation is decision based
        % depending on the order (for efficiency)
        while(~(fleq(alpha, 0.5) || fleq(alpha, 0.1) || fleq(alpha, 0)))
            alpha = alpha - 0.2;
            v02 = v02 + 1;
        end
        
        % 0.1 / 0.5 part approximation
        if ~fleq(alpha, 0)
            G = G * carlapp(fix_s(1/alpha), N, Gi, get_init_gain(b, a, sign(alpha_init)*alpha), retol);
        end
        
        % 0.2 part approximation
        if ~fleq(v02, 0)
            G = G * carlapp(5, N, Gi, get_init_gain(b, a, sign(alpha_init)*0.2), retol)^v02;
        end
        
        % Identification error
        err = 0;
        
    else
        
        % Alpha is smaller than allowed resolution
        if alpha < (1/M)

            if round(alpha*M) >= 1/M

               % Alpha can be rounded to minimal resolution
               G = carlapp(M, N, Gi, get_init_gain(b, a, sign(alpha_init)*(1/M)), retol);
               err = abs(round(alpha)-alpha);
               
            else
                
               % Alpha too small, thus negligible
               G = zpk(1);
               err = alpha;
               
            end
            
        else

            % Use decomposition method with emphasis on efficiency
            % Iterate through the fractions
            for P=2:M      
                
                if alpha >= (1/P)

                    % Get approximation
                    G = G * carlapp(P, N, Gi, get_init_gain(b, a, sign(alpha_init)*(1/P)), retol);
                    alpha = alpha - (1/P);
                    
                end
            
            end
            
            % Error
            err = alpha;
            
        end

    end

    % Add integer part, if any
    G = G * Gi^v;
    
    % Final model reduction
    if ~isempty(retol)
        G = minreal(G, retol);
    end
    
    % Identification quality index
    J   = imp_freq_index(G, Gi_init, alpha_init, w);
    
end

% Compute the initial gain for Carlson approximation of the
% fractional power zero-pole pair
function G0 = get_init_gain(b, a, alpha)

    
    wm = 1/sqrt(a*b);                          % Compute center frequency
    G0 = abs((b*wm*1i+1)/(a*wm*1i+1))^alpha;   % Compute magnitude at wm
    
end

