function [wt, ws, wo] = csens(Gct, A, B, wi)
%CSENS Fractional system sensitivity function frequencies.
%
%      Computes the frequencies for the complementary senitivity function
%      and sensitivity function for the fractional control system given
%      by Gct(jw)=G(jw)*C(jw).
%
%      The complementary sensitivity function is given by
%
%              |  C(jw)*G(jw)  |
%      T(jw) = | ------------- | [dB]
%              | 1+C(jw)*G(jw) |
%
%      And the sensitivity function is given by
%
%              |       1       |
%      S(jw) = | ------------- | [dB]
%              | 1+C(jw)*G(jw) |
%
%      Frequencies sought adhere to the rules:
%      For all w >= wt [rad/sec] -> | T(jwt) | = A [dB]
%      For all w <= ws [rad/sec] -> | S(jws) | = B [dB]
%
%      Usage:   [WT, WS, W] = CSENS(GCT, A, B)
%
%      where    WT, WS - frequencies discussed above,
%               W - array with the whole frequency range [rad/s],
%               A<0, B<0 - noise/disturbance attenuation in dB.
%
%               [WT, WS] = CSENS(GCT, A, B, W)
%               
%      See also: fotf, freqresp, bode, logspace

    % Load configuration parameters
    config = fomcon('config');

    if nargin < 3
        error('CSENS:NotEnoughInputArguments', ...
              'Not enough input arguments.');
    end

    % Default values
    minExp = config.Core.Frequency_domain_computations.Min_freq_exp;
    maxExp = config.Core.Frequency_domain_computations.Max_freq_exp;
    numPts = config.Core.Frequency_domain_computations.Num_points;
    
    if nargin > 3
        lspace = wi;
    else
        % Create logspace object
        lspace=logspace(minExp,maxExp,numPts);
        wo = lspace;
    end
    
    % Check input values
    if A > 0
        warning('CSENS:dBPositiveNoiseAttenuation', ...
                'Given noise attenuation is positive: A > 0 [dB]');
    end
    
    if B > 0
        warning('CSENS:dBPositiveSensFun', ...
                'Given sensitivity function is positive: B > 0 [dB]');
    end
    
    % Complementary sensitivity function
    T = feedback(Gct, 1);
    
    % Sensitivity function
    S = inv(1+Gct);
    
    % Use fsolve to obtain solutions for wt and ws
    wt = find_freq(T, A, lspace);
    ws = find_freq(S, B, lspace);

end

% Computes the response in dB
function w = find_freq(G, a, lspace)

    rsp = freqresp(G, lspace);    
    c = abs(abs(20*log(abs(rsp))/log(10)) + a*ones(1,length(lspace)));
    [ig1,i] = min(c);
    w = lspace(i);
    
end

