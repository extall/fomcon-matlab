function G1 = trunc(G, na_acc, a_acc)
%TRUNC Truncate the exponents and coefficients of a fractional-order transfer function.
%
%      Usage: G1 = trunc(G, na_acc, a_acc)
%      where       G - initial system
%                  na_acc - exponent accuracy (i.e. 1e-2 means 2 decimal
%                  places)
%                  a_acc - coefficient accuracy

    [a,na,b,nb,ioDel] = fotfparam(G);
    
    if  na_acc ~= 0
        na = fix_s(na/na_acc)*na_acc;
        nb = fix_s(nb/na_acc)*na_acc;
    end
    
    if nargin == 3 && a_acc ~= 0
        % Also truncate coefficients
        a = fix_s(a/a_acc)*a_acc;
        b = fix_s(b/a_acc)*a_acc;
    end
    
    G1 = simple(fotf(a,na,b,nb,ioDel));

end

