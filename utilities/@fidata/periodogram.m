function varargout = periodogram(id, opt)
%PERIODOGRAM Returns periodogram of the input and output signals
%   USAGE: [pxx,f] = periodogram(id, opt)
%     where id - identification dataset,
%           opt - options, can be 'u', 'y', or 'both' (by default, shows
%                 periodograms for both input and output, i.e., 'both')
%

% Check input arguments
if nargin < 1
    error('Not enough input arguments.');
end

% Make sure default value is used for option if not supplied
if nargin < 2 || (~strcmpi(opt,'u') && ~strcmpi(opt, 'y') ...
        && ~strcmpi(opt,'both'))
    opt = 'both';
end

% Get the values (TODO: make MIMO safe)
u = id.u;
y = id.y;
fs = 1/id.dt;

% Do the comps
[up, f, w] = compperiodogram(u,fs);
[yp, f, w] = compperiodogram(y,fs);

% Based on the comps, assign output and draw the plots
if (nargout < 1)
    h = figure;
end


switch lower(opt)
    case 'both'
        
        outp = [colv(up) colv(yp)];
        h1 = subplot(211); plot(w, up); title('Input periodogram');
        ylabel('Power/Frequency [dB/rad/s]'); grid;
        h2 = subplot(212); plot(w, yp); title('Output periodogram');
        xlabel('Frequency [rad/s]'); ylabel('Power/Frequency [dB/rad/s]');
        grid;
        linkaxes([h1,h2], 'x');
        set(h1,'xscale','log');
        set(h2,'xscale','log');
      
    case 'y'
        
        outp = colv(y);
        plot(w, yp); title('Input periodogram');
        set(gca,'xscale','log');
        ylabel('Power/Frequency [dB/rad/s]');
        xlabel('Frequency [rad/s]');
        grid;
        
    case 'u'
        
        outp = colv(u);
        plot(w, up); title('Input periodogram');
        set(gca,'xscale','log');
        ylabel('Power/Frequency [dB/rad/s]');
        xlabel('Frequency [rad/s]');
        grid;

end

% Assign outputs
if nargout == 1
    varargout{1} = outp;
elseif nargout == 2
    varargout{2} = f;
end

end

function [dbhz, freq, w] = compperiodogram(x,Fs)
    
    % Use FFT to compute the periodogram
    N = length(x);
    xdft = fft(x);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    
    % Final result in dB/Hz
    dbhz = 10*log10(psdx);
    
    % Get freq, convert Hz to rad/s
    freq = 0:Fs/length(x):Fs/2;
    w = hz2rads(freq);
end