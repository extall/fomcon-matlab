function [polystr] = fpoly2str(a, na, baseVar, mulSymbol)
%FPOLY2STR     Converts a polynomial array set to a string
%
% Usage:
%
%   POLYSTR = FPOLY2STR(A, NA) where A is term coefficients vector
%             and NA is the corresponding term exponents vector
%
%   POLYSTR = FPOLY2STR(A, NA, BASEVAR) uses BASEVAR as base variable
% 
%   POLYSTR = FPOLY2STR(A, NA, BASEVAR, MULSYMBOL) in addition uses the 
%             MULSYMBOL as the multiplication symbol (it is omitted by
%             default, i.e., MULSYMBOL = "")
%
%   Note: Coefficient and exponent arrays are sorted according to the
%         values in the exponent array in descending order prior to
%         generate a polynomial string
%
%   Examples:
%
%        fpoly2str([1 -2 3], [0.5 1 0])      returns '-2s+s^0.5+3'
%        fpoly2str([5 1 -9], [1.5 7 3], 'z') returns 'z^7-9z^3+5z^1.5'   
%
%   See also: str2fpoly
    
	% Load FOMCON configuration (remove this if using elsewhere!)-|
	config = fomcon('config');                                   %|
	numSigDig = config.Core.General.Model_significant_digits;    %|
	
    % Get base variable
    if nargin >= 3
        varn = baseVar;
    else
        varn = 's';
    end

    % Get multiplication symbol
    if nargin >= 4
        muls = mulSymbol;
    else
        muls = '';
    end
    
    % Sort input arrays
    polyArray = [a' na'];
    polyArray = sortrows(polyArray,-2);
    
    % Form a string
    polystr = '';
    
    for n=1:size(polyArray,1)
       % Term coefficient
       wasCoeff = 0;  % Determine if there was a coefficient to append the multiplication symbol
       if ~fleq(polyArray(n,1),1) && ~fleq(polyArray(n,1),-1) || fleq(polyArray(n,2), 0)
           wasCoeff = 1;
           if polyArray(n,1) > 0 && n ~= 1
               polystr=strcat(polystr,'+',num2str(polyArray(n,1),numSigDig));
           else
               polystr=strcat(polystr,num2str(polyArray(n,1),numSigDig));
           end
       elseif fleq(polyArray(n,1),-1)
           polystr=strcat(polystr,'-');
       elseif fleq(polyArray(n,1),1) && n ~= 1
           if ~fleq(polyArray(n,1), 1) && ~fleq(polyArray(n,1), -1)
               wasCoeff = 1;
           end
           polystr=strcat(polystr,'+');
       end
  
       % Term exponent
       if fleq(polyArray(n,2),1)
           % Multiplication symbol
           if wasCoeff
                polystr = strcat(polystr, muls);
           end
           polystr=strcat(polystr,varn);
       elseif polyArray(n,2) ~= 0
           if wasCoeff
                polystr = strcat(polystr, muls);
           end
           polystr=strcat(polystr, varn, '^', num2str(polyArray(n,2),numSigDig));
       end
       
    end
    
end
