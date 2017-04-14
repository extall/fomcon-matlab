function [r, w] = freqresp(G, w)
%FREQRESP  Frequency response of fractional-order transfer functions.
%
%   Usage: R = freqresp(G, W)
%   where
%          G - fractional-order transfer function,
%          w - frequencies, at which the frequency response must be
%          computed.

config = fomcon('config');

if nargin < 1
    error('FREQRESP:NotEnoughInputArguments', ...
        'Not enough input arguments.');
end

minExp = config.Core.Frequency_domain.Default_min_freq_exp;
maxExp = config.Core.Frequency_domain.Default_max_freq_exp;
numPts = config.Core.Frequency_domain.Default_num_points;

if nargin < 2
    w = logspace(minExp,maxExp,numPts);
end

a = G.a;
na= G.na;
b = G.b;
nb= G.nb;
j=sqrt(-1);

for i=1:length(w)
    P=b*((j*w(i)).^nb.');
    Q=a*((j*w(i)).^na.');
    r(i)=P/Q;
end

% Delay
if G.ioDelay > 0
    
    for i=1:length(w)
        r(i) = r(i) * exp(-j*w(i)*G.ioDelay);
    end
    
end

end
