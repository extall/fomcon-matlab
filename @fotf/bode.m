function varargout = bode(G,w)
% BODE Bode frequency response of fractional dynamic systems.
%
% Usage:    [MAG, PH] = BODE(G, W)
%           [MAG, PH, W] = BODE(G)
%           H = BODE(G)
%           H = BODE(G,W)
%           BODE(G)
%           BODE(G, W)
%
% where
%           MAG - magnitude
%           PH - phase
%           W  - frequency range [rad/s]
%           H  - FRD object
%
%           If no output argument is supplied, BODE will create a plot
%           with the frequency response.
%
%           See also: freqresp, frd


if nargin < 1
    error('BODE:NotEnoughInputArguments', ...
        'Not enough input arguments.');
end

% Load configuration parameters
config = fomcon('config');

minExp = config.Core.Frequency_domain.Default_min_freq_exp;
maxExp = config.Core.Frequency_domain.Default_max_freq_exp;
numPts = config.Core.Frequency_domain.Default_num_points;

if nargin==1
    w=logspace(minExp,maxExp,numPts);
end

H1=frd(freqresp(G, w),w);

if nargout==0
    bode(H1);
elseif nargout==2
    % Get magnitude and phase
    [varargout{1}, varargout{2}] = bode(H1, w);
elseif nargout==3
    % Get magnitude, phase and w
    [varargout{1}, varargout{2}, varargout{3}] = bode(H1);
else
    varargout{1} = H1;
end

end