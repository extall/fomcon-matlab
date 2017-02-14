function varargout = validate(id, G)
%VALIDATE Validate fractional model time-domain identification results
%
% Usage:  ERR = VALIDATE(ID, G|FSIM)
%
%         where ERR is the optional output argument which contains the
%                   error vector y(t) - y_id(t) at time points
%                   specified in ID. If omitted, plots validation
%                   results
%
%               ID is a FIDATA structure with validation data
%
%               G is the model to validate or a FSPARAM structure if
%               validation is done on the basis of approximation
%
% See also: fidata, fid

    % Initialize output arguments
    varargout = {};

    % Number of arguments
    if nargin < 2
        error('VALIDATE:NotEnoughInputArguments', ...
            'Not enough input arguments');
    end
    
    % Determine the type of simulation (GL or approximated classic LTI)
    if ~isa(G, 'fsparam')
        % Simulate system using GL based solver
        y = lsim(G, id.u, id.t);
    else
        % Use the approximation
        Z = oustapp(G.plant, G.w(1), G.w(2), G.N, G.approx);
        y = lsim(Z, id.u, id.t);
    end
    

    % Calculate error
    err = y - id.y;
    
    % Fitness measure
    fitness = 100*(1-(norm(err))/(norm(id.y-mean(id.y))));
    
    % Check no. of output arguments
    if nargout == 0
        
        % Plot results
        h = figure();

        subplot(3,1,1);
        plot(id.t, id.y, id.t, y, '--', 'Linewidth', 2);
        ylabel('System output y(t)');
        legend('Initial data','Identified model','Location','Best');
        grid;

        subplot(3,1,2);
        plot(id.t, id.u, 'Linewidth', 2);
        ylabel('System input u(t)');
        grid;
        
        subplot(3,1,3);
        plot(id.t, err, 'Color', 'r', 'Linewidth', 2);
        xlabel('Time [s]');
        ylabel('Output error');
        grid;

        set(h, 'NumberTitle', 'off');
        set(h, 'Name', 'Time-domain validation results');
        title(['Error norm: ' num2str(sum(err.^2)) '; Fit: ' sprintf('%.2f',fitness) '%']);
        
    elseif nargout > 0
        
        varargout{1} = err;
        
    end     
    
end

