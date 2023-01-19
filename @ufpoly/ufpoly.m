classdef ufpoly
    %UFPOLY Fractional polynomial with uncertainty in coeffs and exponents
    %
    % Calling syntax: p = ufpoly(a, na)
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
    % p = ufpoly('
    %
    %  
    %   
    
    properties
        Property1
    end
    
    methods
        function obj = upoly(inputArg1,inputArg2)
            %UPOLY Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

