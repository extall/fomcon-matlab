function [polystr] = poly2str(a, na, baseVar)
%POLY2STR     Converts a polynomial array set to a string
%
% Usage:
%
%   POLYSTR = POLY2STR(A, NA) where A is term coefficients vector
%             and NA is the corresponding term exponents vector
%
%   POLYSTR = POLY2STR(A, NA, BASEVAR) uses BASEVAR as base variable
%
%   Note: Coefficient and exponent arrays are sorted according to the
%         values in the exponent array in descending order prior to
%         generate a polynomial string
%
%   Examples:
%
%        poly2str([1 -2 3], [0.5 1 0])      returns '-2s+s^0.5+3'
%        poly2str([5 1 -9], [1.5 7 3], 'z') returns 'z^7-9z^3+5z^1.5'   
%
%   See also: str2poly
    
	% Load FOMCON configuration (remove this if using elsewhere!)-|
	config = fomcon('config');                                   %|
	numSigDig = config.Core.General.Model_significant_digits;    %|
	
    % Get base variable
    if nargin ==3
        varn = baseVar;
    else
        varn = 's';
    end
    
    % Sort input arrays
    polyArray = [a' na'];
    polyArray = sortrows(polyArray,-2);
    
    % Form a string
    polystr = '';
    
    for n=1:size(polyArray,1)
       % Term coefficient
       if ~fleq(polyArray(n,1),1) && ~fleq(polyArray(n,1),-1) || fleq(polyArray(n,2), 0)
           if polyArray(n,1) > 0 && n ~= 1
               polystr=strcat(polystr,'+',num2str(polyArray(n,1),numSigDig));
           else
               polystr=strcat(polystr,num2str(polyArray(n,1),numSigDig));
           end
       elseif fleq(polyArray(n,1),-1)
           polystr=strcat(polystr,'-');
       elseif fleq(polyArray(n,1),1) && n ~= 1
           polystr=strcat(polystr,'+');
       end
       
       % Term exponent
       if fleq(polyArray(n,2),1)
           polystr=strcat(polystr,varn);
       elseif polyArray(n,2) ~= 0
           polystr=strcat(polystr, varn, '^', num2str(polyArray(n,2),numSigDig));
       end
       
    end
    
end
