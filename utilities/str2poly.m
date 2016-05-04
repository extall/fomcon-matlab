function [c, f] = str2poly(polystr)
%STR2POLY Converts a string to a polynomial
%
%   [C, F] = STR2POLY(POLYSTR) Converts an input string POLYSTR with 
%   terms of base variable 's' to a set of polynomial term coefficients C
%   and respective exponents F. Fractional exponents are also accepted.
%
%   Note: Coefficient and exponent arrays are sorted according to the
%         values in the exponent array in descending order
%
%   Example:
%
%       [c, f]=str2poly('15s^3+2s+1')     returns c=[15 2 1] and f=[3 1 0]
%
%   See also: poly2str
    
    % Replace brackets and remove spaces
    polystr = strrep(polystr, '{', '(');
    polystr = strrep(polystr, '}', ')');
    polystr = strrep(polystr, '[', '(');
    polystr = strrep(polystr, ']', ')');
    polystr = strrep(polystr, ' ', '');
    
    % Locate missing '*' characters
    mulMissing = regexp(polystr,'(?m)(\-?\+?[\.0-9]+([eE]\-?[0-9]+)?)(s(\^\(?\-?\+?[\.\/0-9]+([eE]\-?[0-9]+)?)?)\)?(?=(\+|\-|$|\(|\)))', 'tokens');
    
    for n=1:length(mulMissing)
       
       term = mulMissing{n};
       
       % Add '*' to term
       polystr = strrep(polystr, [term{1} term{2}], [term{1} '*' term{2}]);
       
    end
    
    % Resolve bracket multiplication
    polystr = strrep(polystr, ')(', ')*(');
    
    % Base variable
    s=fotf('s');
    
    G_poly = fotf(eval(polystr));
    [c,f]=fotfparam(G_poly);
    
end
