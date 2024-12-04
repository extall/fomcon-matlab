function is_stable = mikhailovfoclcs(Controller, Plant, w)
% MIKHAILOVFOCLCS - Analyze the stability of a closed-loop fractional-order system using the Mikhailov criterion.
%
% Syntax:
%   is_stable = mikhailovfoclcs(Controller, Plant, w)
%
% Inputs:
%   Controller - A structure containing numerator and denominator polynomials of the controller:
%   Plant      - A structure containing numerator and denominator polynomials of the plant:
%   w          - (Optional) Frequency range for analysis. Default: -1000:0.01:1000.
%
% Outputs:
%   is_stable  - Complex frequency response values (real and imaginary parts).
%
% Description:
%   This function computes the Mikhailov curve for a fractional-order closed-loop system 
%   defined by the given `Controller` and `Plant`. It normalizes the system using the 
%   highest-order terms and plots the Nyquist diagram for visualization.
%   Stability is assessed based on the trajectory of the Mikhailov curve: 
%   encirclements of the origin indicate instability.
%
% Example:
%   % Define controller and plant structures
%Plant=ufotf('1.3s^.3 + 1.4', '15s^1.6 + 2.5s^.3+ 1.5',.1)
%Controller = ufotf('2s^.3 +.5','s^.3');
%mikhailovfoclcs(Controller, Plant)
% If frequency range is not provided, use default frequency range
 if nargin < 3
     w = -1000:0.01:1000;
 end
%High order of the denaminator of Controller
HOC = Controller.a.na(1,1)
%The coefficient of high order of the denaminator of Controller
CHOC = Controller.a.a(1,1)
%High order of the denaminator of Plant
HOP = Plant.a.na(1,1)
%The coefficient of high order of the denaminator of Plant
CHOP = Plant.a.a(1,1)
syms('s');
Ca = ufpoly2str(Controller.a,'*');
Ca = eval(Ca)
Cb = ufpoly2str(Controller.b,'*');
Cb = eval(Cb)
Pa = ufpoly2str(Plant.a,'*');
Pa = eval(Pa)
Pb = ufpoly2str(Plant.b,'*');
Pb = eval(Pb)
Cl = Ca*Pa+Cb*Pb*exp(-Plant.ioDelay(:,1)*s)

highcoe = CHOC*CHOP;
malpha= HOC+HOP;

dcp = Cl/(highcoe*(s+1)^malpha)

 fr = matlabFunction(dcp);

% Compute the response
r = fr(sqrt(-1)*w); rr = real(r); ri = imag(r);
disp(size(r))
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
% 
% 
