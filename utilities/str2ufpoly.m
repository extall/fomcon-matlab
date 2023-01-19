function [a, na, symb] = str2ufpoly(poly)
%STR2UFPOLY String with fractional polynomial with uncert. to vectors
% Calling sequence:
%  [A, NA, SYMB] = STR2UFPOLY(POLY)
%
% The function is based on a parser for the advanced syntax to read
% polynomials with uncertainty intervals defined in square brackets.
%
% Inputs:
% POLY - the string containing the polynomial to process. 
% The intervals themselves are specified with square brackets in the string
% and the lower and upper bounds can be delimited with either , (comma) or
% ; (semicolon). The symbol in the polynomial is inferred automatically and
% returned along with the other parameters.
%
% Outputs:
% A - matrix with as many rows as there are coefficients and two columns
% for the lower and upper bounds of the uncertainty interval.
% NA - the same sized matrix that holds the bounds for all exponents.
% SYMB - the variable symbol. NB! Do not use common variable symbols that
% denote specific items such as e, i, j.
%
% Example:
%
%  [a, na, s] = str2ufpoly('-3s^[1.3,1.5] + [12; 15] s^1.3e-9 - [9, 11]')

if nargin < 1 || ~isa(poly, 'char')
    error('Please provide a string.');
end

% Remove all spaces and perform other replacements
poly = regexprep(poly,'[\s]','');
poly = strrep(poly, '{', '(');
poly = strrep(poly, '}', ')');

% Lowercase only.
poly = lower(poly);

% Automatically infer the variable symbol. Use of "e" is forbidden.
myloc = isstrprop(poly, 'alpha');
symbols = poly(myloc);
symbols(symbols=='e') = [];
if any(myloc) && isempty(symbols)
    error('Please do not use the symbol ''e'' to denote the variable.');
end

% Check if the there is a variable symbol (free term could have been
% supplied)...
if isempty(symbols)
    % Use default symbol 's'
    symb = 's';
else
    symb = symbols(1);
end

% ...and if the symbol is unique
symbols(symbols==symb) = [];
if ~isempty(symbols)
    error('The inferred polynomial variable is not unique.')
end

% Insert multiplication symbols as necessary
poly = placemuls(poly, symb);

% Split the terms by signs. Two special cases need
% to be accounted for: 
% 1. First sign.
% 2. Scientific notation as in 4e+10 or 1e-5. 
% Both are remedied with temporary replacements of these tokens.

% Instead of parsing, tokenize everything that can prevent from isolating
% terms of the polynomial
poly = exp2token(poly);
poly = intervalsigns2token(poly); % Signs inside intervals

if strcmp(poly(1), '+') || strcmp(poly(1), '-')
    poly = [sign2token(poly(1)) poly(2:end)];
end

% Break the polynomial into separate terms by sign
signlocs = find((poly == '+') | (poly == '-'));
signlocs = [1 signlocs length(poly)+1];

terms = {};
for k=1:length(signlocs)-1
    terms{k} = poly(signlocs(k):signlocs(k+1)-1);
end

numTerms = length(terms);

% Initialize the parser
a = zeros(numTerms, 2);
na = zeros(numTerms, 2);

% Iterate over the terms and extract the values depending on what they are
for k=1:numTerms
    term = terms{k};

    % Bring back the correct expressions
    term = token2exp(term);
    term = token2intervalsigns(term);
    term = token2sign(term);

    [a1, na1] = parseTerm(term, symb);
    a(k,:) = a1; na(k,:) = na1;

end


end

% Parse a single term and get the information
function [a, na] = parseTerm(term, symb)

    a = zeros(1, 2); na = zeros(1, 2);

    if isFreeTerm(term, symb)
        % Parse this as a free term, meaning na remains [0, 0].
        a = parseInterval(term);
    else
        % Break into two parts according to mul. sign and parse the two
        % possible number locations. If no mul. sign, then coeff. is
        % simply 1.
        term_parts = explode_str_in_two(term, '*');
        
        if isempty(term_parts)
            a = [1 1];
            na_part = term;
        else
            a = parseInterval(term_parts{1});
            na_part = term_parts{2};
        end

        % Now deal with na_part, which is related to the exponent of the
        % specific variable. Here, instead of a mul sign, we check for
        % exponentiation sign ^. If no such sign, na is simply [1 1].
        exp_parts = explode_str_in_two(na_part, '^');
        if isempty(exp_parts)
            na = [1 1];
        else
            na = parseInterval(exp_parts{2});
        end

    end

end

% Parse the content of [a, b] or [a; b]
function bounds = parseInterval(str)

    if haveIntervals(str)

        % Store the sign. Need to multiply the resulting
        % array by sign later, otherwise sign information is lost.
        % TODO: This actually depends on whether it is agreed by default
        %  that sign for uncertainty intervals is always '+' and it is
        %  the content of the intervals that determines the sign.
        tsign = 1;
        if strcmp(str(1), '-')
            tsign = -1;
        end

        % First attempt
        interval_parts = explode_str_in_two(str, ',');

        % Second attempt
        % TODO: Make more elegant, can use regexp
        if isempty(interval_parts)
            interval_parts = explode_str_in_two(str, ';');
        end

        if isempty(interval_parts)
            error(['Invalid interval specified: ' str]);
        end

        % Proceed with extraction of data
        first = interval_parts{1};
        second = interval_parts{2};

        % Remove sign if necessary
        if strcmp(first(1), '+') || strcmp(first(1), '-')
            first = first(2:end);
        end

        % Remove [ and ]
        first = first(2:end);
        second = second(1:end-1);

        % And just convert to numbers
        lb = str2double(first);
        ub = str2double(second);

        % Multiply by sign
        lb = tsign*lb;
        ub = tsign*ub;

        if lb > ub
            warning(['Interval bounds in ' str ' must be lb <= ub, reversing order.']);
            rev = lb; lb = ub; ub = rev;
        end

    else
        % Very simple case.
        lb = str2double(str);
        ub = str2double(str);
    end

    bounds = [lb, ub];
    
end

% This function is unnecessary in newer releases of MATLAB,
% but for compatibility it is implemented here as well.
function parts = explode_str_in_two(str, delimiter)
    loc = find(str == delimiter);
    parts = {};
    if ~isempty(loc)
        parts{1} = str(1:loc-1);
        parts{2} = str(loc+1:end);
    end
end

function res = isFreeTerm(str, symb)
    res = ~any(find(str == symb));
end

function res = haveIntervals(str)
    res = any(find(str == '['));
end

% Helping functions for properly breaking the polynomial into terms.
% TODO: Using these functions may be less effective than implementing a
%  fully fledged polynomial parser. Thus, please consider implementing a
%  parser if this becomes an issue in the sequel.

function str1 = exp2token(str)
    str1 = str;
    str1 = strrep(str1, 'e-', 'EXPMINUS');
    str1 = strrep(str1, 'e+', 'EXPPLUS');
end


function str1 = token2exp(str)
    str1 = str;
    str1 = strrep(str1, 'EXPMINUS', 'e-');
    str1 = strrep(str1, 'EXPPLUS', 'e+');
end

function str1 = intervalsigns2token(str)
    str1 = str;
    str1 = strrep(str1, '[-', '[MINUS');
    str1 = strrep(str1, '[+', '[PLUS');
    str1 = strrep(str1, ';-', 'SEMICOLMINUS');
    str1 = strrep(str1, ',-', 'COMMAMINUS');
    str1 = strrep(str1, ';+', 'SEMICOLPLUS');
    str1 = strrep(str1, ',+', 'COMMAPLUS');
end

function str1 = token2intervalsigns(str)
    str1 = str;
    str1 = strrep(str1, '[MINUS', '[-');
    str1 = strrep(str1, '[PLUS', '[+');
    str1 = strrep(str1, 'SEMICOLMINUS', ';-');
    str1 = strrep(str1, 'COMMAMINUS', ',-');
    str1 = strrep(str1, 'SEMICOLPLUS', ';+');
    str1 = strrep(str1, 'COMMAPLUS', ',+');
end

function str1 = token2sign(str)
    str1 = str;
    str1 = strrep(str1, 'MINUS', '-');
    str1 = strrep(str1, 'PLUS', '+');
end

function str1 = sign2token(str)
    str1 = str;
    str1 = strrep(str1, '-', 'MINUS');
    str1 = strrep(str1, '+', 'PLUS');
end

% Placing multiplication signs where necessary
function newstr = placemuls(oldstr, symb)

    newstr = oldstr;


    % Find locations of symbols. If any is not preceded by *, insert
    locs = find(newstr == symb);

    if ~isempty(locs)
        % First symbol means the string begins with 's ...' and there is no
        % multiplication there necessary
        if locs(1) == 1
            locs = locs(2:end);
        end
    
        added = 0;
        for k=1:length(locs)
            if ~strcmp(newstr(locs(k)-1), '*')
                % Need to insert *
                newstr = [newstr(1:locs(k)-1+added) '*' newstr(locs(k)+added:end)];
                added = added + 1;
            end
        end
        
        % Resolve bracket multiplication
        newstr = strrep(newstr, ')(', ')*(');  
    end

end

