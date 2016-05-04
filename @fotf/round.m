function G1 = round(G, na_acc, a_acc)
%ROUND Rounds the exponents and coefficients of a FOTF system
%
%      Usage: G1 = round(G, na_acc, a_acc)
%      where       G - initial system
%                  na_acc - exponent accuracy (i.e. 1e-2 means 2 decimal
%                  places)
%                  a_acc - coefficient accuracy

    [a,na,b,nb,ioDel] = fotfparam(G);
    
    if  na_acc ~= 0
        na = round(na/na_acc)*na_acc;
        nb = round(nb/na_acc)*na_acc;
    end
    
    if nargin == 3 && a_acc ~= 0
        % Also round coefficients
        a = round(a/a_acc)*a_acc;
        b = round(b/a_acc)*a_acc;
    end
    
    G1 = simple(fotf(a,na,b,nb,ioDel));

end

