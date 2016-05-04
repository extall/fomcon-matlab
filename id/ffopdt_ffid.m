function [K,L,T,alpha,steps] = ffopdt_ffid(w,k,Kdc,L,H)
%FFPDT_FID Frequency-domain identification of FO-FOPDT model
%
%          [K,L,T,ALPHA] = FFPDT_FID(W,K,KDC,L,H) identifies a
%          fratcional-order FOPDT model of the form
%
%                                    Kdc        -L*s
%                        G(s) = -------------- e
%                                        alpha
%                               1 + T * s
%
%          to be used in auto-tuning. The user must supply a set of
%          identified points {(w1,k1), (w2,k2), ...} where each entry 
%          corresponds to a single (frequency, gain) point on the
%          frequency response plot. The parameters T and alpha are then
%          determined. The delay value must also be supplied but is not
%          used in the computations at this point, so the same value will 
%          be returned by the function. At least two points are needed.
%          H must be a set of strictly decreasing step values for
%          computation, e.g. H=[0.25 0.1 0.01] (default). ALPHA range is
%          permanently set to [0, 2] and cannot be changed.

    % Check number of parameters
    if nargin < 4
        error('FFOPDT:NotEnoughInputArguments', ...
              'Not enough input arguments.');
    end
    
    if nargin < 5
        H = [0.25 0.1 0.01];
    end
    
    % Make sure step size is positive
    H = abs(H);
    
    % DEBUG
%     aa = 0.5:0.01:1;
%     JJ=zeros(length(aa));
%     for pa=1:length(aa)
%        JJ(pa) = compute_cost(w,k,Kdc,aa(pa)); 
%     end
%     h1=figure;
%     plot(aa,JJ);
%     title('Alpha vs J');
%     disp('Hit any key to cont w algorithm');
%     pause;
%     close(h1);
    
    % *** BEGIN ALGORITHM ***
    
    stepcount = 0;
    
    % Initial alpha
    alpha = 1;
    
    % Set up step counter
    Nh        = length(H);
    kh        = 1;
    
    % *** First step ***
    
    h1        = H(kh); kh = kh + 1;
    J_initial = [compute_cost(w,k,Kdc,alpha-h1);
                 compute_cost(w,k,Kdc,alpha);
                 compute_cost(w,k,Kdc,alpha+h1)];
             
    % Minimum value
    [minJ, indJ] = min(J_initial);

    % *** Second step ***
    
    % Determine search direction, positive by default
    h_sign = 1;
    if (indJ == 1), h_sign = -1; end
    
    % *** Step 2.1: determine general point ***
    
    % Number of steps
    k_N = floor((1-h1)/h1);
    
    for m=0:k_N
       alpha_test      = alpha+(m*h_sign*h1);
       alpha_next_test = alpha+((m+1)*h_sign*h1);
       alpha_found     = alpha_test;
       if compute_cost(w,k,Kdc,alpha_test) < ...
          compute_cost(w,k,Kdc,alpha_next_test)
          break;
       end
       stepcount = stepcount + 1;
    end
    
    % Set new alpha value
    alpha = alpha_found;

    % *** Step 2.2: search in the neighborhood of the point ***
    while (kh <= Nh)

        % Determine first value
        alpha = (alpha - h_sign * H(kh-1)) - h_sign * H(kh);
        
        % Number of steps
        k_N = floor(2*(H(kh-1))/H(kh));
        
        for m=0:k_N
        
            alpha_test      = alpha+(m*h_sign*H(kh));
            alpha_next_test = alpha+((m+1)*h_sign*H(kh));
            alpha_found     = alpha_test;
            
            if compute_cost(w,k,Kdc,alpha_test) < ...
               compute_cost(w,k,Kdc,alpha_next_test)
            break;
            end
            
            stepcount = stepcount + 1;
        end
        
        alpha = alpha_found;
        
        % Increase step counter
        kh=kh+1;
        
    end
    
    % *** Third step ***
    
    % All done, return Kdc and compute T
    K = Kdc;
    T = mean(compute_time_constant(w,k,Kdc,alpha));
    steps = stepcount;
end


% Compute the time constant given the identification point (K*,w*),
% the static gain of the process, and the fractional order thereof
function T = compute_time_constant(w,A,K,alpha)

    % Parameterization in constant values of K and alpha
    Tp = @(w,A) (sqrt(K^2-sin(alpha*pi/2)^2.*A.^2)-A.*cos(alpha*pi/2))./...
                (A.*w.^alpha);
    T = Tp(w,A);
    
end

% % Approximate way to compute the time constant
% function T = approx_compute_time_constant(ws,Ks,K,alpha)
%     T = K / (ws^alpha * Ks);
% end

% Compute the cost function J(alpha)
function J = compute_cost(w,k,Kdc,alpha)
    
    %                        ^
    % Compute all entries in T
    That = compute_time_constant(w,k,Kdc,alpha);
    T    = mean(That);

    % Magnitude response
    ffopdt_abs = @(w) sqrt(Kdc^2 ./ (1+2*cos(alpha*pi/2).*w.^alpha.*T+w.^(2*alpha).*T.^2));
    
    pp = ffopdt_abs(w);
    
    % The cost is given by a sum of squares
    J = sum((pp-k).^2);

end