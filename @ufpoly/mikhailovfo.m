function is_stable = mikhailovfo(P, w)
%safop Check the stability of a FO control system

% If frequency range is not provided, use default frequency range
if nargin < 2
    w = -1000:0.01:1000;
end

% Closed loop system
%CL = feedback(C*G, 1);
%G = 
% The highest order in s
%[malpha, ind] = max(P.na(:,1))
malpha = max(P.na(:,1));
highcoe = P.a(1,1);
% Extract the characteristic polynomial
%[cp,~, ~] = fotf2expr(P);
cps = ufpoly2str(P,'*');
% Perform the division
syms('s');
cp = eval(cps)
dcp = cp/(highcoe*(s+1)^malpha)

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


