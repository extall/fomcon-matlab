function [Kp, Ki, Kd, lam, mu, Z] = fpid_optimize_(fsim, fopt, G, varargin)
%FPID_OPTIMIZE_ Obtain optimal values for FOPID control loops
%   This is a wrapper function for the fpid_optimize() function. Instead of
%   outputting the log of the underlying optimization algorithm to the
%   command line, this function creates a graphical user interface which
%   displays progress information related to the optimization progress.
%   The calling sequence is identical to fpid_optimize(), so please see the
%   documentation for that function to see the argument list.
%
% See also: fpid_optimize, fsparam, fpopt

% Get the relevant optimization options out of the structure
opt_type   = fopt.type;
opt_p      = fopt.p;
opt_metric = fopt.metric;
opt_alg    = fopt.alg;
opt_op     = fopt.optop;

% Create the additional variables to be passed to the output function
addVars = struct;
addVars.p    = opt_p;
addVars.type = opt_type;
addVars.alg  = opt_alg;

% Designate the output function
opt_op.OutputFcn = @(x, optimValues, state) ...
    fpid_optimize_OutputFcn(x, optimValues, state, addVars);
opt_op.Display = '';

% Re-specify the optimization options
new_fopt       = fopt;
new_fopt.optop = opt_op;

% Close the window if already opened
fpid_optimize_pui_close();

% Specify the algorithm
myAlg = fpid_optimize_capitalize_words(opt_alg);

% SQP get special treatment here
if strcmpi(myAlg, 'sqp'), myAlg = upper(myAlg); end

% Formulate the final string
alg_str = [myAlg ' / ' upper(opt_metric)];

d = guidata(fpid_optimize_pui_upd());

% Set up the labels
switch lower(opt_alg)
    
    case 'interior-point'
        set(d.lblOptimalityType, 'String', 'First order optimality');
        
    case 'sqp'
        set(d.lblOptimalityType, 'String', 'First order optimality');
        
    case 'active-set'
        set(d.lblOptimalityType, 'String', 'First order optimality');
        set(d.lblFeasibility, 'String', 'Max constraint');
        set(d.lblStepNorm, 'String', 'Directional derivative');
        
    case 'nelder-mead'
        
        set(d.txtStepNorm, 'Enable', 'Off');
        set(d.txtFeasibility, 'Enable', 'Off');
        
end

% Set up the optimized parameters
switch lower(opt_type)
    
    case 'n'
        
        
    case 'e'
        set(d.txtLambda, 'Enable', 'Off');
        set(d.txtMu, 'Enable', 'Off');
        
    case 'g'
        set(d.txtKp, 'Enable', 'Off');
        set(d.txtKi, 'Enable', 'Off');
        set(d.txtKd, 'Enable', 'Off');
        
end

% Clear parameter values
set(d.txtAlgorithm, 'String', alg_str);
set(d.txtIterationNum, 'String', '');
set(d.txtPerformance, 'String', '');
set(d.txtSimulationNum, 'String', '');
set(d.txtOptimalityType, 'String', '');

if strcmpi(opt_alg, 'nelder-mead')
    set(d.txtStepNorm, 'String', 'N/A');
    set(d.txtFeasibility, 'String', 'N/A');
else
    set(d.txtStepNorm, 'String', '');
    set(d.txtFeasibility, 'String', '');   
end

set(d.txtKp, 'String', num2str(opt_p(1)));
set(d.txtKi, 'String', num2str(opt_p(2)));
set(d.txtLambda, 'String', num2str(opt_p(4)));
set(d.txtKd, 'String', num2str(opt_p(3)));
set(d.txtMu, 'String', num2str(opt_p(5)));

drawnow();

% Run the algorithm
[Kp, Ki, Kd, lam, mu, Z] = fpid_optimize(fsim, new_fopt, G, varargin{:});

% Save the final performance log to GUI as structure
h = fpid_optimize_pui_upd();
stopSignal = getappdata(h,'StopSignal');

% Inform user that identification has stopped
if stopSignal<1
    d = guidata(h);
    set(d.txtIdentStatus, 'String', 'Optimization terminated.');
    fid_pui_deactivate_actions();
end

% If the operation was aborted, assign previous parameter values
if (stopSignal == 2)
    Kp  = opt_p(1);
    Ki  = opt_p(2);
    Kd  = opt_p(3);
    lam = opt_p(4);
    mu  = opt_p(5);
end

it = getappdata(h, 'Iteration');
perfLog = getappdata(h, 'PerformanceLog');
paramLog = getappdata(h, 'ParameterLog');
if it>=1
    perfLog = perfLog(1:it);
    paramLog = paramLog(1:it,:);
end
elapsedTime = getappdata(h, 'Toc');

savedPerfData = struct;
savedPerfData.metric = upper(opt_metric);
savedPerfData.iterations = 1:it;
savedPerfData.performance = perfLog;
savedPerfData.parameters = paramLog;
savedPerfData.elapsed_time = elapsedTime;

setappdata(h, 'savedPerfData', savedPerfData);

end

% Function that creates and maintains the GUI for the pUI
function h = fpid_optimize_pui_upd()
    % Locate the pUI in the list of currently open GUI windows
    h = findall(0, 'Tag', 'figFpidOptimizePui');
    % If there is no such window, create it, and store
    % an empty performance graph log in there: initially 50 values
    
    if isempty(h)
        h = fpid_optimize_pui();
        setappdata(h,'PerformanceLog', zeros(50,1));
        setappdata(h,'ParameterLog', zeros(50,5));
        setappdata(h,'Iteration',0);
        curTic = tic; setappdata(h, 'Tic', curTic);
        curToc = toc(curTic); setappdata(h,'Toc', curToc);

    end
end

% Disable the buttons once identification completes
function fid_pui_deactivate_actions()
h = fpid_optimize_pui_upd();
d = guidata(h);
set(d.btnStop, 'Enable', 'off');
set(d.btnAbort, 'Enable', 'off');
end

function fpid_optimize_pui_close()
    h = findall(0, 'Tag', 'figFpidOptimizePui');
    if ~isempty(h)
        close(h);
    end
end

% OutputFcn for the optimization algorithm
function stop = fpid_optimize_OutputFcn(x, optimValues, state, addVars)

% Get additional variables
opt_type = addVars.type;
opt_alg  = addVars.alg;

% Stop is FALSE by default
stop = false;

% Collect iteration data
iteration = optimValues.iteration;
perfIndex = optimValues.fval;
simCount  = optimValues.funccount;

% Create or update the process UI figure
h = fpid_optimize_pui_upd();

% Fetch optimized parameters based on optimization type
x_kp     = addVars.p(1); x_ki     = addVars.p(2);
x_lambda = addVars.p(4); x_kd     = addVars.p(3);
x_mu     = addVars.p(5);

switch lower(opt_type)
    
    case 'n'
        
        x_kp     = x(1); x_ki     = x(2);
        x_lambda = x(4); x_kd     = x(3);
        x_mu     = x(5);
        
    case 'e'
        
         x_kp     = x(1);
         x_ki     = x(2);
         x_kd     = x(3);
        
    case 'g'
        
        x_lambda = x(1);
        x_mu     = x(2);
        
end

if (iteration>0)
    
    % Output the necessary data
    d = guidata(h);

    % Update the controller parameters
    set(d.txtKp, 'String', num2str(x_kp));
    set(d.txtKi, 'String', num2str(x_ki));
    set(d.txtLambda, 'String', num2str(x_lambda));
    set(d.txtKd, 'String', num2str(x_kd));
    set(d.txtMu, 'String', num2str(x_mu));
    
    set(d.txtIterationNum, 'String', num2str(iteration));
    set(d.txtPerformance, 'String', num2str(perfIndex));
    set(d.txtSimulationNum, 'String', num2str(simCount)); 
    
    % Fetch optimization data based on the algorithm and update the fields
    switch lower(opt_alg)
        
        case {'interior-point', 'sqp'}
            set(d.txtOptimalityType, 'String', ...
                num2str(optimValues.firstorderopt));
            set(d.txtStepNorm, 'String', ...
                num2str(optimValues.stepsize));
            set(d.txtFeasibility, 'String', ...
                num2str(optimValues.constrviolation));
        
        case 'active-set'
            set(d.txtOptimalityType, 'String', ...
                num2str(optimValues.firstorderopt));
            set(d.txtStepNorm, 'String', ...
                num2str(optimValues.directionalderivative));
            set(d.txtFeasibility, 'String', ...
                num2str(optimValues.constrviolation));
            
        case 'nelder-mead'
            
            set(d.txtOptimalityType, 'String', optimValues.procedure);
            
    end
    
    % Fetch the performance log
    perfLog = getappdata(h,'PerformanceLog');
    if iteration > size(perfLog,1)
        perfLog = [perfLog; zeros(50, 1)];
    end
    perfLog(iteration) = perfIndex;
    
    % Fetch the parameter log
    paramLog = getappdata(h,'ParameterLog');
    if iteration > size(paramLog,1)
        paramLog = [paramLog; zeros(50, 5)];
    end
    paramLog(iteration,1:5) = [x_kp x_ki x_kd x_lambda x_mu];
    
    % Save the performance graph related information
    setappdata(h,'PerformanceLog', perfLog);
    setappdata(h,'ParameterLog', paramLog);
    setappdata(h,'Iteration', iteration);
    curTic = getappdata(h, 'Tic');
    setappdata(h, 'Toc', toc(curTic));
    
    % Update figure window
    drawnow();
    
end

% Check for stop condition
if (getappdata(h,'StopSignal')>=1)
    stop = true;
end

end

function cap_str = fpid_optimize_capitalize_words(str)
 cap_str = regexprep(str, '(\<[a-z](?!\s))', '${upper($0)}');
end