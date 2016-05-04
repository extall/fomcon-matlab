function [a,na,b,nb,T] = fotfparam(G)
% FOTFPARAM  get FOTF object parameters
%
% Usage; [A,NA,B,NB,T] = FOTFPARAM(fotf_object)
%
%        Returns FOTF object parameters where
%        A, NA - pole polynomial coefficients and exponents,
%        B, NB - zero polynomial coefficients and exponents,
%        T     - ioDelay [sec]
%
    if nargout == 2
        a = G.b;
        na = G.nb;
    else
        a  = G.a;
        na = G.na;
        b  = G.b;
        nb = G.nb;
        T = G.ioDelay;
    end
end
