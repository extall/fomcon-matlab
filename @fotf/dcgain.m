function K = dcgain(G)
%DCGAIN Static gain of the fractional transfer function
%
%   K = dcgain(G) computes the steady-state gain of the fractional-order
%                 system G.

    [a,na,b,nb] = fotfparam(G);
    a(~fleq(na,0))=0; b(~fleq(nb,0))=0; % Use only static coefficients
    K = polyval(b,0)/polyval(a,0);

end

