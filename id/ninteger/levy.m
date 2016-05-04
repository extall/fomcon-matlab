function [G, J, handle] = levy (w, gain, phase, Q, n, m)

% [G, J, handle] = levy (w, gain, phase, Q, n, m)
% This function returns a model for a plant
% with a known frequency behaviour.
% The model is
% b(1)*s^(m*Q) + b(2)*s^((m-1)Q) + ... + b(end-2)*s^(2*Q) + b(end-1)*s^Q + b(end)
% -------------------------------------------------------------------------------
% a(1)*s^(n*Q) + a(2)*s^((n-1)Q) + ... + a(end-2)*s^(2*Q) + a(end-1)*s^Q + a(end)
% The function receives the gain in dB and the phase in degrees
% for frequencies, in rad/s, specified in w.
% Parameters Q, n and m must be supplid;
% a structure is returned with fields num and den,
% containing vectors b and a respectively.
% J is the quadratic error per sampling frequency.
% If the last output argument exists, a handle for a Bode plot will be returned.
% Duarte Valério 2003

if size(w, 1) < size(w, 2)
    w = w'; % now w is a column vector for sure
end
if size(gain, 1) < size(gain, 2)
    gain = gain'; % now w is a column vector for sure
end
if size(phase, 1) < size(phase, 2)
    phase = phase'; % now w is a column vector for sure
end

R = real(10.^(gain/20) .* exp(j * deg2rad(phase)));
I = imag(10.^(gain/20) .* exp(j * deg2rad(phase)));

A = [];
for l = 0:m
    for c = 0:m
        A(l+1,c+1) = sum(-real((j*w).^(l*Q)).*real((j*w).^(c*Q))...
            -imag((j*w).^(l*Q)).*imag((j*w).^(c*Q)));
    end
end

B = [];
for l = 0:m
    for c = 1:n
        B(l+1,c) = sum(real((j*w).^(l*Q)).*real((j*w).^(c*Q)).*R...
            +imag((j*w).^(l*Q)).*real((j*w).^(c*Q)).*I...
            -real((j*w).^(l*Q)).*imag((j*w).^(c*Q)).*I...
            +imag((j*w).^(l*Q)).*imag((j*w).^(c*Q)).*R);
    end
end

C = [];
for l = 1:n
    for c = 0:m
        C(l,c+1) = sum( (imag((j*w).^(l*Q)).*I - real((j*w).^(l*Q)).*R) .* real((j*w).^(c*Q))...
            +(-real((j*w).^(l*Q)).*I - imag((j*w).^(l*Q)).*R) .* imag((j*w).^(c*Q)) );
    end
end

D = [];
for l = 1:n
    for c = 1:n
        D(l,c) = sum(real((j*w).^(l*Q)).*real((j*w).^(c*Q)).*(R.^2+I.^2)...
            +imag((j*w).^(l*Q)).*imag((j*w).^(c*Q)).*(R.^2+I.^2));
    end
end

e = [];
for l = 0:m
    e(l+1,1) = sum(-real((j*w).^(l*Q)).*R -imag((j*w).^(l*Q)).*I);
end

f = [];
for l = 1:n
    f(l,1) = sum(-real((j*w).^(l*Q)).*(R.^2+I.^2));
end

coefficients = [A B; C D] \ [e; f];
G.num = coefficients(m+1:-1:1);
G.den = [coefficients(end:-1:m+2); 1];

% index J
if nargout > 1
    tempNum = 0;
    for c = 1:length(G.num)
        tempNum = tempNum + G.num(c) * (j*w).^((m-c+1)*Q);
    end
    tempDen = 0;
    for c = 1:length(G.den)
        tempDen = tempDen + G.den(c) * (j*w).^((n-c+1)*Q);
    end
    J = sum((abs(10.^(gain/20).*exp(j*deg2rad(phase)) - tempNum./tempDen)).^2) / length(w);
end

% Bode plots
if nargout == 3
    wNew = logspace(floor(log10(w(1))), ceil(log10(w(end))), max(50, length(w)));
    tempNum = 0;
    for c = 1:length(G.num)
        tempNum = tempNum + G.num(c) * (j*wNew).^((m-c+1)*Q);
    end
    tempDen = 0;
    for c = 1:length(G.den)
        tempDen = tempDen + G.den(c) * (j*wNew).^((n-c+1)*Q);
    end
    gainNew = 20 * log10(abs(tempNum ./ tempDen));
    phaseNew = rad2deg(unwrap(angle(tempNum ./ tempDen)));
    handle = figure;
    subplot(2,1,1)
    semilogx(w,gain,'.', wNew,gainNew)
    subplot(2,1,2)
    semilogx(w,phase,'.', wNew,phaseNew)
end