classdef ufotf
    %UFOTF Fractional-order continuous transfer function with uncertainties
    
    properties(SetAccess = protected)
        b        % Fractional zeros polynomial with uncertainty intervals
        a        % Fractional poles polynomial -"-
        ioDelay  % Delay that can be represented by an uncertainty interval
    end
    
    methods
        function obj = ufotf(b, a, ioDelay)
            %UFOTF Create a UFOTF object
            
            % Check input arguments
            if nargin < 2
                error('To create a UFOTF object, please provide UFPOLYs (or strings) B and A.');
            end

            if ~isa(b, 'ufpoly') && ~isa(b, 'char')
                error('B must either be a UFPOLY object or string.');
            end

            if ~isa(a, 'ufpoly') && ~isa(a, 'char')
                error('A must either be a UFPOLY object or string.');
            end

            if nargin < 3
                ioDelay = [];
            end

            if nargin == 3 && (all(size(ioDelay) ~= [1,1]) && ...
                    all(size(ioDelay) ~= [1, 2]))
                error('The delay parameter, if specified, must be a scalar, or a 1x2 vector');
            end

            if nargin == 3 && min(ioDelay) < 0
                error('The ioDelay parameter must be a positive scalar or contain an interval in [0, inf)');
            end

            % Assign the parameters
            if isa(b, 'ufpoly')
                obj.b = b;
            elseif isa(b, 'char')
                obj.b = ufpoly(b);
            end

            if isa(a, 'ufpoly')
                obj.a = a;
            elseif isa(a, 'char')
                obj.a = ufpoly(a);
            end

            if nargin == 3
                if isscalar(ioDelay)
                    obj.ioDelay = [ioDelay ioDelay];
                else
                    obj.ioDelay = ioDelay;
                end
            else
                obj.ioDelay = [];
            end


        end
        
    end
end

