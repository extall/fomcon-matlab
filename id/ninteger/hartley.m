function G = hartley (w, gain, phase, Q, n, type)

% G = hartley (w, gain, phase, n, ordmax, type)
% This function returns a model for a plant
% with a known frequency behaviour.
% The model is either
% G(1)*s^(n*Q) + G(2)*s^((n-1)Q) + ... + G(end-2)*s^(2*Q) + G(end-1)*s^Q + G(end)
% or else the inverse of the above.
% The function receives the gain in dB and the phase in degrees
% for frequencies, in rad/s, specified in w.
% Parameter Q must be supplid;
% if n (the desired number of zeros and poles of the model)
% is not supplied, length(w)-1 is assumed instead;
% type may be 'num' or 'den' for the two types of models stated above:
% it may be ommited, the default being 'den'.
% Duarte Valério 2003

if nargin < 5
    n = length(w) - 1;
end
if nargin < 6
    type = 'den';
end
if size(w, 1) < size(w, 2)
    w = w'; % now w is a column vector for sure
end
if size(gain, 1) < size(gain, 2)
    gain = gain'; % now w is a column vector for sure
end
if size(phase, 1) < size(phase, 2)
    phase = phase'; % now w is a column vector for sure
end

if strcmp(type, 'num')
    g = 10.^(gain/20) .* exp(j * deg2rad(phase));
else
    g = 1 ./ (10.^(gain/20) .* exp(j * deg2rad(phase)));
end
W = [];
for i = 0:n
    W = [W, (j*w).^(i*Q)];
end
if n==length(w)
    k = W\g;
else
    k = pinv(W)*g;
end
G = k(end:-1:1);