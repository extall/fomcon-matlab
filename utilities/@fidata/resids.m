function varargout = resids(id, G, opts)
%VALIDATE Validate fractional model time-domain identification results
%
% Usage:  [ERR, STATS] = RESIDS(ID, G|FSIM, OPTS)
%
%         where ERR is the optional output argument which contains the
%                   absolute error vector y(t) - y_id(t) at time points
%                   specified in ID.
%               STATS is a structure with the following entries:
%                     .MaxError -- maximum absolute value of the simulation
%                                  error,
%                     .Mean -- mean value of residuals,
%                     .ResidNorm -- residual norm: sum(err^2)
%                     .MSE -- mean squared error: 1/N sum(err^2)
%                     .RErr -- vector with autocorrelation values of
%                              residuals for lags tau up until tau_max
%
%               If output arguments are omitted, plots validation results
%
%         Inputs:
%               ID is a FIDATA structure with validation data,
%               G is the model to generate statistics for, or a FSPARAM
%               structure if residual analysis is done on the basis of
%               approximation,
%               OPTS is a structure with the following entries (optional)
%                    .Conf -- measure of confidence in (0.0, 1.0]. Defaults
%                             to 0.95 (i.e., 95% confidence),
%                    .MaxTau -- maximum number of lags to consider for the
%                               autocorrelation of residuals test, defaults
%                               to 50
%
% See also: fidata, fid, validate

    % Initialize output arguments
    varargout = {};

    % Number of arguments
    if nargin < 2
        error('VALIDATE:NotEnoughInputArguments', ...
              'Not enough input arguments');
    end
    
    % Default values
    maxtau = 50;
    conf = 0.95;
    
    % Check optional arguments
    if nargin == 3
         
        if ~isa(opts, 'struct')
            error('RESIDS:OptionsNotAStructure', ...
                  'The OPTS argument must be a structure');
        end
        
        % Confidence
        if cfieldexists(opts, 'Conf')
            conf = opts.Conf;
        end
        
        % Maximum lags
        if cfieldexists(opts, 'MaxTau')
            maxtau = opts.MaxTau;
        end
        
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

    % maxtau must be less than the amount of simulated points
    if length(y)<=maxtau
        maxtau = length(y)-1;
    end
    
    % Calculate simulation error
    err = y - id.y;
    
    % Confidence interval
    intc = quantile(c_p(conf))/sqrt(length(err));
    
    % Compute necessary parameters
    resnorm = sum(err.^2);              % Residual norm
    mse     = 1/length(err)*resnorm;    % Mean squared error
    errmean = mean(err);                % Mean value
    [maxerror, ind] = max(abs(err));    % Maximum error
    
    % Residuals at lags
    rerr = [];
    for k=0:maxtau
        rerr(end+1)=R_e(err,k);
    end
    
    % V1: why?
    %remax = max(abs(rerr));
    %rerr=rerr/remax;
    
    % By first lag
    rerr = rerr/rerr(1);
    
    % Check no. of output arguments
    if nargout == 0
        
        % Plot results
        h = figure();
        
        subplot(2,1,1);
        plot(id.t, err, 'Color', 'r', 'Linewidth', 2);
        hold on;
        plot([id.t(1) id.t(end)], [errmean errmean], '--b', 'LineWidth', 2);
        stem(id.t(ind), err(ind), 'k', 'LineWidth', 2);
        xlabel('Time [s]');
        ylabel('Output error');
        title(['Mean squared error: ' num2str(mse) '; Max abs error: ' num2str(maxerror)]);
        grid;
        
        subplot(2,1,2);
        stem(1:maxtau,rerr(2:end), 'Color', 'b', 'Linewidth', 2);
        hold on;
        plot([1 maxtau], [intc intc], ':r', 'Linewidth', 2);
        plot([1 maxtau], -[intc intc], ':r', 'Linewidth', 2);
        xlabel('Lags [Samples]');
        xlim([1 maxtau]);
        title(['Autocorrelation of residuals (with P=' num2str(conf) ' confidence)']);

        set(h, 'NumberTitle', 'off');
        set(h, 'Name', 'Residual analysis');
    
    end
    
    if nargout > 0    
        varargout{1} = err;
    end     
    
    if nargout > 1
        stats = struct;
        stats.MaxError = maxerror;
        stats.Mean = errmean;
        stats.ResidNorm = resnorm;
        stats.MSE = mse;
        stats.RErr = rerr;
        varargout{2} = stats;
    end
    
end

% Residual covariance
function Re = R_e(err,tau)
    Re = 1/(length(err)-tau)*sum(err(1:end-tau).*err(1+tau:end));
end

% Argument for quantile function: assuming normal distribution
function cp = c_p(p)
    cp = 1-0.5*(1-p);
end

