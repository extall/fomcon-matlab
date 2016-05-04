function n = get10exp(x)
%GET10EXP Get n from 10^n == x
%  Usage: n=get10exp(x)
    
    n = fix_s(log(x)/log(10));

end

