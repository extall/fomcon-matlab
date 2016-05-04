function J = imp_freq_index(H, G, alpha, w, doPlot)
% IMP_FREQ_INDEX Implicit fractional transfer function approximation index calculation

    % Load configuration parameters
    config = fomcon('config');
    numPts = config.Core.Frequency_domain_computations.Num_points;

    % Check input argument number
    if nargin < 5, doPlot = false; end
    
    if nargin < 4
        error('IMPFREQINDEX:NotEnoughInputArguments', ...
              'Not enough input arguments.');
    end

    % Frequency range
    wb = w(1);
    wh = w(2);
    
    % Take N log-spaced points
    w = logspace(log10(wb), log10(wh), numPts);
    
    % Original response
    r_orig = squeeze(freqresp(G,w)).^alpha;
    
    % Approximation
    r = squeeze(freqresp(H,w));
    
    % Performance index
    J = sum(abs(r_orig - r).^2)/length(r);
    
    % Also plot BODE diagram if requested
    if doPlot
        bode(frd(r_orig, w),w);
        hold on;
        bode(frd(r, w),w);
        legend('Original response', 'Obtained approximation');
    end

end