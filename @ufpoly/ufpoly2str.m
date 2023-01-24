function str = ufpoly2str(p, mul)
%UFPOLY2STR Convert UFPOLY to string

if nargin < 2
    mul = '';
end

str = ufpoly2str(p.a, p.na, p.symb, mul);

end

