function [yr,ur,w] = get_fft_fidata(id1)
%GET_FFT_FIDATA Get the FFT for input/output data for FIDATA object (test)

% Fetch the data
y1 = id1.y;
u1 = id1.u;
dt = id1.dt;

% Get both signals' FFT
yr = do_fft(y1,dt);
[ur,w] = do_fft(u1,dt);

end

function [r,w] = do_fft(y,dt)
    Fs = 1/dt;
    L = length(y);
    NFFT = 2^nextpow2(L);
    Y = fft(y,NFFT)/L;
    w = hz2rads(Fs/2*linspace(0,1,NFFT/2+1)); % Convert to rad/s
    r = Y(1:NFFT/2+1);
end
