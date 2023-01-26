function [polystr, uintvl] = ufpoly2str(a, na, baseVar, mulSymbol)
%UFPOLY2STR Converts a fractional polynomial with uncert. into string
%
% Usage:
%
%   POLYSTR = UFPOLY2STR(A, NA) where A is term coefficients matrix
%             and NA is the corresponding term exponents matrix, with
%             uncertainties, so for both the dimensions are N x 2
%
%   POLYSTR = UFPOLY2STR(A, NA, BASEVAR) uses BASEVAR as base variable
% 
%   POLYSTR = UFPOLY2STR(A, NA, BASEVAR, MULSYMBOL) in addition uses the 
%             MULSYMBOL as the multiplication symbol (it is omitted by
%             default, i.e., MULSYMBOL = "")
%
%   See also: fpoly2str
    
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
    
    % Form a string
    polystr = '';
    
    % Count uncertainty intervals. Will be useful in the display function.
    uintvl = 0;
    
    for n=1:size(a,1)  % Same size, no matter which one to use for iterations
       % Term coefficient
       wasCoeff = 0;  % Determine if there was a coefficient to append the multiplication symbol
       coeffNominalVal = fleq(a(n,1), a(n,2));
       expNominalVal = fleq(na(n,1), na(n,2));
       expUnity = fleq(na(n,1), 1);
       expZero = fleq(na(n,1),0) && fleq(na(n,2),0);
       
       % Coefficient
       if coeffNominalVal
          coeff = a(n,1);
          
          % Is it unity?
          add_sign = '';
          add_coeff = '';
          if ~fleq(abs(coeff), 1)
              add_coeff = num2str(coeff, numSigDig);
              wasCoeff = 1;
          end
          
          % Account for the case of unity as free term
          if fleq(abs(coeff), 1) && expZero
             add_coeff = '1';
          end

          % Account for a negative free term
          if fleq(abs(coeff), 1) && expZero && coeff < 0
             add_sign = '-';
          end
          
          if n == 1 && fleq(abs(coeff), 1) && coeff < 0
              add_sign = '-';
              wasCoeff = 1;
          end
          
          if n ~= 1
              add_sign = signstr(coeff);
          end
          
          polystr = [polystr add_sign add_coeff];
              
       else
           % Just typeset interval adding "+" unless first coeff
           if n ~= 1
               polystr = [polystr '+'];
           end
           polystr = [polystr typeset_interval(a(n, :), numSigDig)];
           wasCoeff = 1;
           uintvl = uintvl + 1;
       end
       
       % Exponent
       if expZero
           % Do nothing
       else
           addmul = '';
           if wasCoeff
               addmul = muls;
           end
           if expNominalVal
               
               if expUnity
                   polystr = [polystr addmul varn];
               else
                   polystr = [polystr addmul varn '^{' ...
                   num2str(na(n,1), numSigDig) '}'];
               end
               
           else
               polystr = [polystr addmul varn '^{' ...
                   typeset_interval(na(n,:), numSigDig) '}'];   
               uintvl = uintvl + 1;
           end
       end
       
    end
    
end

function str = typeset_interval(a, numSigDig)
    str = ['[' num2str(a(1), numSigDig) ...
        ', ' num2str(a(2), numSigDig) ']'];
end

function str = signstr(a)
    str = '+';
    if a < 0, str = '-'; end
end
