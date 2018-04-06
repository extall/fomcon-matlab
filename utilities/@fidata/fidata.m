classdef fidata
    %FIDATA  Time-domain identification parameters for FOTF systems.
    %
    % Usage:  P = FIDATA(Y, U, T, OPTIONS)
    %
    % where
    %         Y - observed system output
    %         U - observed system input
    %         T - observations time vector
    %         OPTIONS - various options
    %
    % Note:   If parameters are not equally sized, then resizing or
    %         trimming will occur in the following order: Y, U, T.
    %
    % The resulting object contains correctly sized identification
    % parameters and additionally sampling interval DT, based on the
    % provided time vector (or sample time). All parameters can be accessed
    % via usual dot notation, e.g. to get the output vector type
    %
    %         y = fidata_object.y
    % 
    % Note that parameter names are case sensitive!
    %
    % See also: fid, fidata/plot, fidata/validate
    
    properties
        y       % Observed system output
        u       % Observed system input
        t       % Observations time vector
        dt      % Sample time
        tstmp   % Timestamp showing date/time of dataset creation
    end
    
    methods
        
        % Initializer
        function p = fidata(y, u, t, varargin)
           
            % Check number of input parameters
            if nargin < 3
                error('FIDATA:NotEnoughInputArguments', 'Not enough input arguments.');
            end
                
            % Check vectors - swap to row vectors if needed
            
            % Y
            if size(y,2) ~= 1 && size(y,1) == 1
               
                % Transpose
                y = y';
                
            elseif size(y,1) > 1 && size(y,2) > 1
                
                error('FIDATA:YNotVector', 'Output dataset Y is not a vector!');
                
            end
            
            % U
            if size(u,2) ~= 1 && size(u,1) == 1
               
                % Transpose
                u = u';
                
            elseif size(u,1) > 1 && size(u,2) > 1
                
                error('FIDATA:UNotVector', 'Output dataset U is not a vector!');
                
            end
            
            % T
            if size(t,2) ~= 1 && size(t,1) == 1
               
                % Transpose
                t = t';
                
            elseif size(t,1) > 1 && size(t,2) > 1
                
                error('FIDATA:TNotVector', 'Output dataset T is not a vector!');
                
            end
            
            % Check sizing
            yLen = length(y);
            uLen = length(u);
            tLen = length(t);
            
            % Check data type (must be double)
            if ~isa(y,'double')
                warning('FIDATA:YNotDouble', ...
                        'The Y vector is not of type ''double''. Converting automatically.');
                y = double(y);
            end
            
            if ~isa(u,'double')
                warning('FIDATA:UNotDouble', ...
                        'The U vector is not of type ''double''. Converting automatically.');
                u = double(u);
            end
            
            if ~isa(t,'double')
                warning('FIDATA:YNotDouble', ...
                        'The T vector is not of type ''double''. Converting automatically.');
                t = double(t);
            end
            
            % Determine vector t type
            if tLen == 1
                dt = t;
            else
                % Check whether t is regularly spaced
                if ~fleq(min(diff(t)), max(diff(t)))
                   warning('The time vector is not regularly spaced. Problems may occur during identification.'); 
                end
                
                % Get time step
                dt = t(2) - t(1);
            end

            % Check sizes: Y against U
            if yLen > uLen
                
                % Pad U with last value and issue a warning
                u(end:yLen)=u(end);
                warning('FIDATA:UShorterThanY','Dataset U is shorter than Y, padding U with end value');
            
            elseif yLen < uLen
               
                % Trim U
                u = u(1:yLen);
                warning('FIDATA:ULongerThanY','Dataset U is longer than Y, trimming U');
                
            end
            
            % Refresh length
            uLen = length(u);
            
            % T is sample time
            if tLen == 1
                t = dt*(0:uLen-1)';
                tLen = length(t);
            end
            
            % Check sizes: U against T
            if uLen > tLen
                
                % Resize T
                t(end+1:uLen,:) = t(end,:)+dt*(1:uLen-tLen);
                warning('FIDATA:TShorterThanU','Dataset T is shorter than U, resizing T');
            
            elseif uLen < tLen
                
                % Trim T
                t = t(1:uLen);
                warning('FIDATA:TLongerThanU','Dataset T is longer than U, trimming T');
                
            end
            
            % Set parameters
            p.y = y;
            p.u = u;
            p.t = t;
            p.dt = dt;
            
            % Timestamp
            p.tstmp = datestr(now);
            
        end
        
    end
    
end

