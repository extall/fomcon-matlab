function varargout = sim_regularize(varargin)
%SIM_REGULARIZE Regularize the results of variable-step model simulation
%
% Usage: [Y1, Y2, ..., T] = SIM_REGULARIZE(Y1, Y2, ..., T, DTMIN)
%
% where  Y1, Y2, ... - time domain responses to regularize,
%        T           - non-regularly spaced vector with time samples [s],
%        DTMIN       - minimum allowed step size [s].
%
% The function will compute the time step sizes at all samples and remove
% all Y1, Y2, ... entries at times such that T<DTMIN. A new regularly
% spaced time vector will then be created and Y1, Y2 will be populated with
% values such that Y1(T1:T2) = Y(T1), where T2-T1 is the index of the
% corresponding step size in computation.
%
% Please note that interpolation is not used in this case. The user can
% interpolate the data per the objective using e.g. interp1() function.
%
% See also: fpid_optimize_sim

    % Check minimum number of input arguments
    if nargin < 3
        error('SIM_REGULARIZE:NotEnoughInputArguments', 'Not enough input arguments.');
    end
    
    % Time vector and parameters
    t     = varargin{end-1};
    dtmin = varargin{end};
    
    numY = nargin - 2;                  % Input argument number for system response
    y = cell2mat(varargin(1:numY));     % Fetch initial response
    
    dt = t(2:end) - t(1:end-1);      	% Get time step vector
    
    if max(dt <= dtmin)                 % Determine minimal time step
        dt = dtmin;
    else
        dt = min(dt);
    end

    t_regular = 0:dt:max(t);            % Generate the time vector
    t_ind = fix_s(t/dt);                % Get known key value indicies
    
    % Discard values with zero indicies
    y(t_ind<1,:)=[];
    t_ind(t_ind<1)=[];
    
    t_ind(end+1) = length(t_regular)+1;   % Last element
    
    % Populate the new responses
    y_regular = zeros(length(t_regular),numY);
    for k=1:numY
        for l=1:length(t_ind)-1
            y_regular(t_ind(l):t_ind(l+1)-1,k)=y(l,k);
        end
    end
    
    % Assign outputs
    for k=1:numY, varargout{k}=y_regular(:,k); end    % Responses 
    varargout{numY+1} = t_regular;                    % Time vector
    
end

