function [a,na,b,nb,L,Gid,state]  = fid_(varargin)
%FID_ Fractional-order Transfer Function SISO system identification w/ GUI
%   This is a wrapper function for the fid() function. Instead of
%   outputting the log of the underlying optimization algorithm to the
%   command line, this function creates a graphical user interface which
%   displays progress information related to the identification progress.
%   The calling sequence is identical to fid(), so please see the
%   documentation for that function to see the argument list.
%
% See also: fid, fidata, fotfid, optimset

state='running';

% Build argument list
fidNumArgs = nargin('fid');
in_args = varargin; in_args(end+1:fidNumArgs-1)={[]};

% Check initial system
fs = varargin{1};
G_ini = fotf;
if isa(fs, 'fsparam')
    G_ini = fs.plant;
elseif isa(fs, 'fotf')
    G_ini = fs;
else
    error('FID_:WrongInitialSystem', 'Unexpected initial system!');
end

if nargin < fidNumArgs
    op = optimset;
else
    op = varargin{fidNumArgs};
end

% Designate the output function
op.OutputFcn = @fid_OutputFcn;
op.Display = '';

% Close the window if already opened
fid_pui_close();

% Specify the algorithm used and create the initial window
algName = 'Trust-Region-Reflective';            % Default value
if cfieldexists(op, 'IdentificationAlgorithm')
    switch lower(op.IdentificationAlgorithm)
        case 'trr'
            
        case 'lm'
            algName = 'Levenberg-Marquardt';
    end
end
d = guidata(fid_pui_upd());
set(d.txtAlgorithm, 'String', algName);
set(d.txtIterationNum, 'String', '');
set(d.txtCost, 'String', '');
set(d.txtSimulationNum, 'String', '');
set(d.txtStepSize, 'String', '');
set(d.txtFoOptimality, 'String', '');
drawnow();

% Run the algorithm
[a,na,b,nb,L,Gid] = fid(in_args{1:fidNumArgs-1}, op);

% Save the final performance log to GUI as structure
h = fid_pui_upd();
stopSignal = getappdata(h,'StopSignal');

% Inform user that identification has stopped
if stopSignal<1
    d = guidata(h);
    set(d.txtIdentStatus, 'String', 'Identification terminated.');
    fid_pui_deactivate_actions();
    state = 'stopped';
end

% If the operation was aborted, assign previous parameter values
if (stopSignal == 2)
    [a,na,b,nb,L] = fotfparam(G_ini);
    Gid = G_ini;
    state='aborted';
end

it = getappdata(h, 'Iteration');
perfLog = getappdata(h, 'PerformanceLog');
if it>=1
    perfLog = perfLog(1:it);
end
elapsedTime = getappdata(h, 'Toc');

savedPerfData = struct;
savedPerfData.iterations = 1:it;
savedPerfData.performance = perfLog;
savedPerfData.elapsed_time = elapsedTime;

setappdata(h, 'savedPerfData', savedPerfData);

end

% Function that creates and maintains the GUI for the pUI
function h = fid_pui_upd()
    % Locate the pUI in the list of currently open GUI windows
    h = findall(0, 'Tag', 'figFidPui');
    % If there is no such window, create it, and store
    % an empty performance graph log in there: initially 50 values
    if isempty(h)
        h = fid_pui();
        setappdata(h,'PerformanceLog', zeros(50,1));
        setappdata(h,'Iteration',0);
        curTic = tic; setappdata(h, 'Tic', curTic);
        curToc = toc(curTic); setappdata(h,'Toc', curToc);

    end
end

% Disable the buttons once identification completes
function fid_pui_deactivate_actions()
h = fid_pui_upd();
d = guidata(h);
set(d.btnStop, 'Enable', 'off');
set(d.btnAbort, 'Enable', 'off');
end

function fid_pui_close()
    h = findall(0, 'Tag', 'figFidPui');
    if ~isempty(h)
        close(h);
    end
end

% OutputFcn for the optimization algorithm
function stop = fid_OutputFcn(x, optimValues, state)

% Stop is FALSE by default
stop = false;

% Collect iteration data
iteration = optimValues.iteration;

% Create or update the process UI figure
h = fid_pui_upd();

if (iteration>0)
    cost      = optimValues.resnorm;
    sim_count = optimValues.funccount;
    step_norm = optimValues.stepsize;
    fo_optim  = optimValues.firstorderopt;

    % Output the necessary data
    d = guidata(h);
    set(d.txtIterationNum, 'String', num2str(iteration));
    set(d.txtCost, 'String', sprintf('%0.8e',cost));
    set(d.txtSimulationNum, 'String', num2str(sim_count));
    set(d.txtStepSize, 'String', num2str(step_norm));
    set(d.txtFoOptimality, 'String', num2str(fo_optim));
    
    % Fetch the performance log
    perfLog = getappdata(h,'PerformanceLog');
    if iteration > size(perfLog,1)
        perfLog = [perfLog; zeros(50, 1)];
    end
    perfLog(iteration) = cost;
    
    % Save the performance graph related information
    setappdata(h,'PerformanceLog', perfLog);
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