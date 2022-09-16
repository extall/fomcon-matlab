function dstr = iod2str(L, baseVar, mulSymbol)
%IOD2STR Return a string for a delay term

config = fomcon('config');                                  
numSigDig = config.Core.General.Model_significant_digits;   

% Get base variable
if nargin >= 2
    varn = baseVar;
else
    varn = 's';
end

% Get multiplication symbol
if nargin >= 3
    muls = mulSymbol;
else
    muls = '';
end

if ~isempty(L) && L ~= 0 && L > 0
    dstr = ['exp(-' num2str(L, numSigDig) muls varn ')'];
else
    dstr = '';
end

end

