function r = freqresp(S, w)
%FREQRESP  Frequency response of FO state space models.
%
%   Usage: R = freqresp(S, W)
%   where
%          S - fractional-order state space model,
%          w - frequencies, at which the frequency response is computed.

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

% Get necessary array size
noOuts = size(S.fotf_array,1);
noIns  = size(S.fotf_array,2);

% Create the response array
r = zeros(noOuts,noIns,numel(w));

% Go through all of the FO TFs in the FOSS array
for k=1:noOuts
    for l=1:noIns
        r(k,l,:) = freqresp(getfotf(S,k,l),w);
    end
end

end

