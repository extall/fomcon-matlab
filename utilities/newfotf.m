function  G = newfotf(zeroPoly, polePoly, T)
%NEWFOTF Creates a new FOTF object
%
% Usage:
%      G=NEWFOTF(S1,S2,T) creates a new FOTF objects from provided
%      polynomial strings S1 (zero polinomial) and S2 (pole polynomial)
%      T - optional ioDelay parameter [sec]. Polynomials can also be
%      independently represented by matched numeric vectors
%      [a, na] and [b, nb].
%
%      See also: fotf
    
    if nargin < 2
        error('NEWFOTF:NotEnoughInputArguments', 'Not enough input arguments.');
    end

    if isnumeric(zeroPoly)
        b = zeroPoly(1:end/2);
        nb = zeroPoly(end/2+1:end);
    else
        [b nb] = str2fpoly(zeroPoly);
    end
    
    if isnumeric(polePoly)
        a = polePoly(1:end/2);
        na = polePoly(end/2+1:end);
    else
        [a na] = str2fpoly(polePoly);
    end
    
    if nargin < 3
        T = 0;
    end
    
    G = simple(fotf(a,na,b,nb,T));

end
