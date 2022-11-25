function [c, f] = str2fpoly(polystr, params)
%STR2FPOLY Converts a string to a polynomial
%
%   [C, F] = STR2FPOLY(POLYSTR) Converts an input string POLYSTR with 
%   terms of base variable 's' to a set of polynomial term coefficients C
%   and respective exponents F. Fractional exponents are also accepted.
%   
%   [C, F] = STR2FPOLY(POLYSTR, PARAMS) The optional PARAMS argument makes
%   it possible to use parameters which will simply be replaced to their
%   numerical values. E.g.: params.K=1 will make the STR2POLY function look
%   for the parameter named K in the POLYSTR and if it is found, it will be
%   replaced with its value, i.e., 1.
%
%   Note: Coefficient and exponent arrays are sorted according to the
%         values in the exponent array in descending order
%
%   Example:
%
%       [c, f]=str2fpoly('15s^3+2s+1')     returns c=[15 2 1] and f=[3 1 0]
%
%   See also: fpoly2str
    
    % Replace brackets and remove spaces
    polystr = strrep(polystr, '{', '(');
    polystr = strrep(polystr, '}', ')');
    polystr = strrep(polystr, '[', '(');
    polystr = strrep(polystr, ']', ')');
    polystr = strrep(polystr, ' ', '');
	
	% Search and replace parameters
	if nargin >1 && ~isempty(params)
    
	% Get all parameters
    allModelParams = fieldnames(params);
    for k=1:length(allModelParams)
       polystr = strrep(polystr, allModelParams{k}, ...
             num2str(params.(allModelParams{k})));
    end
	
	end
    
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
