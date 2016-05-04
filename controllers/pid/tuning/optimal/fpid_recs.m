function [p1_lims,p2_lims] = fpid_recs(p1,p2,dp,init,cl_poly_fun,op)
%FPID_RECS Approximate BIBO stability rectangle for FOPID control system
%
%   Usage: [P1_LIMS,P2_LIMS]=FPID_RECS[P1,P2,DP,INIT,CL_POLY_FUN,OP)
%
%   where   P1_LIMS, P2_LIMS - stability limits of first and second coord.
%           
%           P1, P2 - string with first and second coordinates to sweep,
%                    may be one of 'Kp', 'Ki', 'Kd', 'lambda', or 'mu',
%                     NB! All parameter names are case-sensitive!
%           DP     - sweeping step, either scalar or
%                    vector of form [DP1; DP2],
%           INIT   - structure with initial coordinate of the form
%                    INIT.Kp = ..., INIT.Ki = ..., etc.,
%                     NB! Parameter names are defined exactly as above.
%           CL_POLY_FUN - a function handle which must accept five
%                         arguments in the form (Kp,Ki,Kd,lambda,mu) and
%                         return a fotf object having the closed-loop pole
%                         polynomial corresponding to the control system
%                         with the FPID controller
%           OP          - Structure with additional options (optional):
%                         OP.maxPoints = N, where N is the points to check
%                                           from the center point in every
%                                           direction (i.e. from initial
%                                           coordinates), default: 100,
%                         OP.displayPlot = 0 or 1 - draw the approximate 
%                                                   stability region.
%                         OP.progress    = 0 or 1 - display progress of
%                                                   computation
%                         OP.drawUnstable= 0 or 1 - draw unstable points on
%                                                   the stability plane
%                                                   (this only works if
%                                                   displayPlot is
%                                                   activated.)
%
%            The algorithm uses Matignon's stability theorem to determine
%            BIBO stability of closed-loop commensurate-order systems with
%            parallel form FOPID controllers in the loop.
%
%            If systems are not commensurate-order, the stability test
%            my produce unreliable results.
%
% See also: fotf/isstable, fpid_recs_ini

% Check input arguments
if nargin < 5
    error('FPID_RECSTAB:NotEnoguhInputArguments', ...
          'Not enough input arguments.');
end

% Check initial coordinate
allParams = {'Kp', 'Ki', 'Kd', 'lambda', 'mu'};
for n=1:length(allParams)
    if ~cfieldexists(init, allParams{n})
        error('FPID_RECSTAB:MissingInitialParameter', ...
        ['Parameter ' allParams{n} ' not set in initial parameter ' ...
         'structure']);
    end
end

% Check sweeeping parameters
if strcmp(p1,p2)
    error('FPID_RECSTAB:CannotSweepSameParameter', ...
        ['To determine an approximate rectangle it is necessary to ' ...
         'define two parameters to sweep']);
end

% Check whether sweeping parameters are in the initial point
if ~ismember(p1,allParams)
   error('FPID_RECSTAB:BadParameterP1', ...
   ['Parameter ' p1 ' is not a valid FOPID parameter to sweep']);
end

if ~ismember(p2,allParams)
   error('FPID_RECSTAB:BadParameterP2', ...
   ['Parameter ' p2 ' is not a valid FOPID parameter to sweep']);
end

% Check sweeping steps
if length(dp) == 2
    dp1 = dp(1);
    dp2 = dp(2);
elseif length(dp) == 1
    dp1 = dp;
    dp2 = dp;
else
    error('FPID_RECSTAB:WrongSweepingStep', ...
          'Sweeping step must either be scalar or of the form [DP1;DP2].');
end

% *** Check options structure ***

% Max points for sweep
maxPoints = 100;
if nargin == 6 && cfieldexists(op, 'maxPoints')
    maxPoints = op.maxPoints;
end

% Display the plot
displayPlot = 0;
if nargin == 6 && cfieldexists(op, 'displayPlot')
    displayPlot = op.displayPlot;
end

% Display progress
prog = 0;
if nargin == 6 && cfieldexists(op, 'progress')
    prog = op.progress;
end

% Draw unstable points on the plane
dunst = 0;
if nargin == 6 && cfieldexists(op, 'drawUnstable')
    dunst = op.drawUnstable;
end

% *** Main algorithm ***

% Check initial point: it must be stable
if ~isstable(cl_poly_fun(init.Kp, init.Ki, init.Kd, ...
                         init.lambda, init.mu))
    error('FPID_RECSTAB:InitialPointMustBeStable', ...
          'Initial point is unstable.');
end

% Find approximate limits
p1_lims = [];
p2_lims = [];

% Sweep first parameter
p1_lims(1) = do_parameter_sweep(p1, -dp1, maxPoints, init, cl_poly_fun, prog);
p1_lims(2) = do_parameter_sweep(p1, dp1, maxPoints, init, cl_poly_fun, prog);

% Sweep second parameter
p2_lims(1) = do_parameter_sweep(p2, -dp2, maxPoints, init, cl_poly_fun, prog);
p2_lims(2) = do_parameter_sweep(p2, dp2, maxPoints, init, cl_poly_fun, prog);

% Display plot shows approximate stability
% boundaries around the initial point
if displayPlot
    
    % Check the obtained range of parameter values
    x_range = p1_lims(1):dp1:p1_lims(2);
    y_range = p2_lims(1):dp2:p2_lims(2);
    
    % Number of points to check
    numPX = length(x_range);
    numPY = length(y_range);
    
    % Stability matrix
    % s_matrix = zeros(numPY, numPX);
    
    % Check all points and plot stability on the graph
    h = figure;
    
    % Set correct limits
    xlim([p1_lims(1) p1_lims(2)]);
    ylim([p2_lims(1) p2_lims(2)]);
    xlabel(p1);
    ylabel(p2);
    title('Approximate stability plane');
    
    % Plot all points on the same figure
    hold('on');
    
    if prog
       progmsg = 'Visualizing stability plane...';
       progbar = pbar(progmsg);
       maxIter = numPX*numPY;
       curIter = 1;
    end
    
    for k=1:numPY
        for l=1:numPX
            init.(p1) = x_range(l);
            init.(p2) = y_range(k);
            if isstable(cl_poly_fun(init.Kp, init.Ki, init.Kd, ...
                         init.lambda, init.mu))
               plot(x_range(l), y_range(k),'.k');
            else
               if dunst
                  plot(x_range(l), y_range(k),'xr');
               end
            end
            if prog
                % Update the progress bar
                update(progbar, curIter, maxIter, progmsg);
                curIter = curIter + 1;
            end    
        end
    end
    
    % No hold on figure
    hold('off');
    
    % Make the boundaries more defined
    xlims = get(gca, 'xlim');
    ylims = get(gca, 'ylim');
    xlims_new = [xlims(1)-1.5*dp1 xlims(2)+1.5*dp1];
    ylims_new = [ylims(1)-1.5*dp2 ylims(2)+1.5*dp2];
    xlim(xlims_new);
    ylim(ylims_new);
    
    % X/Y color
    sog = 0.3; % Shade of grey
    set(gca, 'xcolor', [sog sog sog]);
    set(gca, 'ycolor', [sog sog sog]);
    
    % Delete progress bar if present
    if prog
        delete(progbar);
    end
    
end
    
end

function [p_lim, stab] = do_parameter_sweep(p, dp, N, params, polyf, pr)
    % Should we display the progress?
    if pr
       progmsg = ['Sweeping ' p ', dp=' num2str(dp) '...'];
       progbar = pbar(progmsg);
    end
    stab = [];
    p_ini = params.(p);
    for k=1:N
        p_lim      = p_ini+k*dp;
        params.(p) = p_lim;
        if ~isstable(polyf(params.Kp, params.Ki, params.Kd, ...
                         params.lambda, params.mu))
            p_lim = p_ini+(k-1)*dp;
            break;
        end
        if pr
            % Update the progress bar
            update(progbar, k, N, progmsg);
        end
    end
    
    % Delete the progress bar
    if pr
        delete(progbar);
    end
end


