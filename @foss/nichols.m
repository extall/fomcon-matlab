function nichols(S,w)
% Nichols frequency response of fractional-order dynamic systems.
%
% See also: nichols

    % Load configuration parameters
    config = fomcon('config');

    minExp = config.Core.Frequency_domain.Default_min_freq_exp;
    maxExp = config.Core.Frequency_domain.Default_max_freq_exp;
    numPts = config.Core.Frequency_domain.Default_num_points;

    if nargin==1
        w=logspace(minExp,maxExp,numPts);
    end
    
    H=bode(S,w);
    nichols(H);
    
end