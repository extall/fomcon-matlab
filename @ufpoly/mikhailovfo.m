function is_stable = mikhailovfo(P, w)
% MIKHAILOVFO - Analyze the stability of a fractional-order system using Mikhailov criterion.
%
% Syntax:
%   is_stable = mikhailovfo(P, w)
%
% Inputs:
%   P - A fractional-order polynomial structure (assumed to contain coefficients and orders).
%   w - (Optional) Frequency range for analysis. Default: -1000:0.01:1000.
%
% Outputs:
%   is_stable - Complex frequency response values (real and imaginary parts).
%
% Description:
%   The function plots the Mikhailov curve and determines stability based
%   on the trajectory of this plot. To assess stability, check the total
%   number of encirclements of the origin by the path. If the number of
%   encirclements equals zero, the polynomial is stable. Otherwise, it is
%   not stable.
% Example:
%   P = ufpoly('s^1.8 + s + s^.65 + 10');
%   mikhailovfo(P)
% If frequency range is not provided, use default frequency range
if nargin < 2
    w = -1000:0.01:1000;
end


malpha = max(P.na(:,1));
highcoe = P.a(1,1);

syms('s');
cps=ufpoly2str(P,'*');
cps=eval(cps)
dcp = cps/(highcoe*(s+1)^malpha);

% Create a function handle for the resulting function
fr = matlabFunction(dcp);

% Compute the response
r = fr(sqrt(-1)*w); rr = real(r); ri = imag(r);
% Plot, if no output argument is provided
if nargout < 1
   
    %plot(rr, ri, 'g', 'LineWidth', 1.5);
    %hold on; plot(0, 0,'r+'); % Origin
end

figure;
f1 = frd(r,w);
nyquistplot(f1);
title('Mikhailov Plot')
hold on
plot(0, 0,'r+')
hold on
plot(-1, 0,'w+')
is_stable = r;


