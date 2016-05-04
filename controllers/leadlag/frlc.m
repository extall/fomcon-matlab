function [r, mag, ph] = frlc(k, x, lam, alpha, w)
%FRLC Fractional Lead-Lag Compensator (FLLC) frequency response
%
% Usage: [R, MAG, PH] = FRLC(K, X, LAM, ALPHA, W)
%
% where   R, MAG, PH  - magnitude and phase responses over W, R - complex
%                       response
%
%         K, X, LAM, ALPHA - FLLC parameters such that
%
%                                              ALPHA
%                              (  LAM*s + 1  )
%                   G(s) = K * (-------------)
%                              ( X*LAM*s + 1 )
%
% Note:   Magnitude and phase returned are those calculated by the BODE
%         function of the Control System toolbox

    % Check number of input arguments
    if nargin < 5
        error('FRLC:NotEnoughInputArguments', 'Not enough input arguments.');
    end

    % Imag unit
    j = sqrt(-1);

    % Calculate response
    r = zeros(length(w),1);
    
    for n=1:length(w)
       num = lam*j*w(n)+1;
       den = x*lam*j*w(n)+1;
       r(n) = k*(num/den)^alpha;
    end
    
    % Only compute magnitude and phase when necessary
    if nargout > 1

        % Get frd object
        h = frd(r, w);

        % Calculate magnitude and phase
        [mag, ph] = bode(h, w);
        
    end

end

