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

            if nargin < 1
                error('Cannot create an empty UFPOLY');
            end

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
                    error('a and na must have exactly two columns (specifying lower and upper bounds, respectively).');
                end

            elseif ~isempty(str)
                % Simply call str2ufpoly
                [a, na] = str2ufpoly(str);
            end
            
            % Remove zero entries, if any
            zero_rows = all(a == 0, 2);
            a(zero_rows) = [];
            na(zero_rows) = [];

            % Sort both matrices according to values in na
            
            % Sort by lower interval value in exponent
            [nal, ind] = sort(na(:,1), 'descend');
            nau = na(ind, 2); na = [nal nau];
            a = a(ind,:);

            % Just set the object properties
            obj.a = a;
            obj.na = na;
            obj.symb = symb;
            
        end
        

    end
end

