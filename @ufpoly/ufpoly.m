classdef ufpoly
    %UFPOLY Fractional polynomial with uncertainty in coeffs and exponents
    %
    % Calling syntax: p = ufpoly(a, na) OR p = ufpoly(str)
    %   where a: coefficient intervals, na: exponent intervals correspon-
    %   ding to a. Sizes must match, rows contain interval boundaries and
    %   the two columns the first must contain the lower boundary, and the
    %   second must contain the upper boundary.
    %
    % It is also possible to pass column vectors to ufpoly, in which case
    % both interval boundaries will be considered equal (the ufpoly then
    % collapses to a fractional polynomial without uncertainties). For that
    % reason, it is also possible to create a ufpoly using the following
    % notation:
    % 
    % p = ufpoly('[-2, -1]s^1.2 + 12s^[0.5, 0.6] - [20, 30]')
    %
    
    properties
        a    % The coefficient bounds for each coefficient
        na   % The exponent bounds for each corresponding exponent
        symb % The variable symbol. 's' by default.
    end
    
    methods
        function obj = ufpoly(varargin)
            %UPOLY Construct an instance of this class
            %   One needs to provide either:
            %   - a, na matrices which hold the coefficient and exponent
            %   bounds
            %   or
            %   - str which holds the string for the ufpoly.
            %   - symb (optional). The variable symbol. Default is 's'.

%             p = inputParser;
% 
%             % The input arguments with default values
%             addOptional(p, 'a', [], @(x) isnumeric(x));
%             addOptional(p, 'na', [], @(x) isnumeric(x));
%             addOptional(p, 'str', '', @ischar);
%             addParameter(p, 'symb', 's', @ischar);
% 
%             parse(p, varargin{:});
% 
%             a = p.Results.a;
%             na = p.Results.na;
%             str = p.Results.str;
%             symb = p.Results.symb;

            symb = 's';
            if isnumeric(varargin{1}) && nargin > 1 && isnumeric(varargin{2})
                a = varargin{1};
                na = varargin{2};
                str = '';
                if nargin > 2
                    symb = varargin{3};
                end
            elseif ischar(varargin{1})
                a = [];
                na = [];
                str = varargin{1};
                if nargin > 1
                    symb = varargin{2};
                end
            end

            if isempty(a) && isempty(na) && isempty(str)
                error('You need to provide at least one of the following: a and na, or str')
            elseif ~isempty(a) && ~isempty(na)
                % Checks on A and NA in terms of shape
                if ~all(size(a)==size(na))
                    error('a and na must have exactly the same dimensions.');
                end
                if size(a, 2) ~= 2 || size(na, 2) ~= 2
                    error('a and na must have exactly two columns.');
                end

            elseif ~isempty(str)
                % Simply call str2ufpoly
                [a, na] = str2ufpoly(str);
            end

            % Sort both matrices according to values in na
            [na, ind] = sort(na, 'descend');
            a = a(ind);

            % Just set the object properties
            obj.a = a;
            obj.na = na;
            obj.symb = symb;
            
        end
        

    end
end

