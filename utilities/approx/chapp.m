function H = chapp(T, alpha, y, wmax)
%CHAPP Charef approximation of an implicit transfer function.
% 
% The function is in form   G(s) = (1 + Ts)^alpha.
%
% Usage: H = CHAPP(T, ALPHA, Y, WMAX)
%
%   where   H - resulting approximation ZPK,
%         
%           T, ALPHA - parameters such that (1+T*s)^ALPHA is
%                         approximated,
%           Y - deviation from the desired response in the
%               frequency range in dB,
%           WMAX - desired system bandwidth in rad/s.

    % Check input arguments
    if nargin < 4
        error('CHAPP:NotEnoughInputArguments', ...
            'Not enough input arguments.');
    end

    % Compute synthesis parameters
    pT = T;
    
    a = 10^(y/(10*(1-alpha)));
    b = 10^(y/(10*alpha));
    ab = 10^(y/(10*alpha*(1-alpha)));
    
    p0 = pT * sqrt(b);
    
    % Iteration count depends on the
    % frequency range of interst
    N = floor(log(wmax/p0)/log(ab)) + 1;
    
    % Create pole/zero arrays
    poles = 1;
    zeros = 1;
    
    s = zpk('s');
    
    % Begin iteration process
    for k=0:N
        poles = poles * (1 + s / (p0 * ab^k));
        if k < N, zeros = zeros * (1 + s / (a * p0 * ab^k)); end
    end
    
    % Obtain approximation
    H = inv(zeros/poles);

end

