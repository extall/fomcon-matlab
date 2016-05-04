function [y, u, t, r] = fpid_optimize_sim(frpid, ltimodel, sopt)
%FPID_OPTIMIZE_SIM Simulate a fractional PID based control system using Simulink.
%
%   Usage: Y = FPID_OPTIMIZE_SIM(FRPID, LTIMODEL, SOPT)
%
%   where
%
%          FRPID    - fractional-order PID controller parameters such that
%                     FRPID = [Kp; Ki; lambda; Kd; mu],
%          LTIMODEL - model of the controlled plant,
%          SOPT   - an options data structure FPID_SIMOPT, containing the following
%                     SOPT.modelname - model name (e.g. 'fpid_optimize_model.mdl'),
%                     SOPT.tmax      - final time value for simulation [s],
%                     SOPT.dtmin     - minimum time step [s],
%                     SOPT.dtmax     - maximum time step [s],
%                     SOPT.ulim      - control law saturation values in form
%                                      [UMIN; UMAX], if empty, the control
%                                      signal is unconstrained
%                     SOPT.r         - reference value (or set value) for
%                                      simulation,
%                     SOPT.w         - Oustaloup approximation frequency
%                                      range [rad/s],
%                     SOPT.N         - Oustaloup filter approximation order,
%                     SOPT.type      - Oustaloup filter approximation type
%                                      ('oust' or 'ref'),
%                     SOPT.cancelzero- Whether to apply zero cancellation
%                                      of the resulting controller to
%                                      ensure that it is proper (0 or 1).
%
%   Outputs:          Y - system output,
%                     U - control signal value,
%                     T - time vector (generally non-regular due to
%                                      variable-size ODE solver employed).
%
%   In order to convert the results to those with a fixed time step (the
%   ODE solver will attempt to find the minimum time step, whereas the user
%   can set the maximum time step), one can use the provided function
%   sim_regularize(). See help for that function for corresponding syntax.
%
% A typical negative unity feedback system is assumed for simulation
%
%         _           ---------------  Satu-   ---------------
%  R    /   \         |             |  ration  |             |    Y
%----->|     | -----> |    GC(s)    |--------->|     G(s)    |----o---->
%    +  \   /         |             |  [UMIN;  |             |    |
%         ^           ---------------   UMAX]  ---------------    | 
%       - |                                                       |
%         ---------------------------------------------------------
%
% NB! the original model 'fpid_optimize_model.mdl' is required for this
% function to work. If you wish to use this model in a different situation,
% please make a copy of it with a different file name rather than editing
% the existing model.
%
% See also: fpid_optimize, sopt, sim_regularize

    % Load configuration parameters
    config = fomcon('config');
	numSigDig = config.Core.General.Model_significant_digits;

    % Check input arguments
    if nargin < 3
        error('FPIDOPTIMIZESIM:NotEnoughInputArguments', ...
              'Not enough input arguments.');
    end
    
    % Fetch fractional PID parameters
    Kp = frpid(1);
    Ki = frpid(2); lam=frpid(3);
    Kd = frpid(4); mu=frpid(5);
    
    % Check input model type
    if ~isa(ltimodel, 'lti')
        error('FPIDOPTIMIZESIM:NotLTIModel', ...
              'Given model is not a LTI model.');
    end
    
    % Check step sizes
    if fleq(sopt.dtmin,sopt.dtmax)
        warning('FPIDOPTIMIZESIM:MinStepEqualToMaxStep', ...
                'Minimum and maximum step sizes are equal.');
        sopt.dtmin = sopt.dtmax - sopt.dtmax/10;
        if sopt.dtmin <= 0
            error('FPIDOPTIMIZESIM:CannotResolveStepSize', ...
                  'Please specify different step sizes for simulation.');
        end
    end
    
    if sopt.dtmin > sopt.dtmax
        warning('FPIDOPTIMIZESIM:MinSizeGreaterMaxSize', ...
                'Minimum step size greater than maximum step size: swapping.');
        m = sopt.dtmin;
        sopt.dtmin = sopt.dtmax;
        sopt.dtmax = m;
    end
    
    % Check control signal limits
    if isempty(sopt.ulim)
        ulim(1) = -Inf;
        ulim(2) = +Inf;
    else
        ulim = sopt.ulim;
    end
    
    % Generate fractional PID controller approximation
    Gc = oustapid(Kp,Ki,lam,Kd,mu,sopt.w(1),sopt.w(2),sopt.N, sopt.type);
    
    % Convert to proper system?
    if sopt.cancelzero
        Gc = toproper(Gc, sopt.w(2));
    end
    
    % Low-pass filter time constant
    tau = 1/sopt.w(2);
    
    % Reference signal
    r = sopt.r;
    
    % What type of reference is it?
    if isa(r, 'refgen')
        u_r  = r.u;
        ts_sim = r.Ts;
    else
        u_r = r;
        ts_sim = sopt.dtmax;
    end
        
    refval = struct;
    refval.signals(1).values = u_r;
    refval.time = [];
    
    % Load model
    modelName = sopt.modelname;
    load_system(modelName);
    
    % Get block names
    blks=find_system(gcs, 'Type', 'block');
    
    % Default system parameters ------------------------------
    % NB! DO NOT change the names of these blocks in the model
    % --------------------------------------------------------
    FRACPID_BLOCKNAME = 'Fractional PID controller';
    PARAMS = {
        {'r'}, {'VariableName'}, {'fpid_optimize_ref'};
        {'r'}, {'SampleTime'}, {num2str(ts_sim, numSigDig)};
        {strcat(FRACPID_BLOCKNAME,'/','FPID LTI')}, {'sys'}, {'fpid_optimize_fpid'};
        {strcat(FRACPID_BLOCKNAME,'/','Lowpass filter 1')}, {'numerator', 'denominator'}, {'1', strcat('[',num2str(tau, numSigDig),' 1]')};
        {'Saturation'}, {'upperlimit', 'lowerlimit'}, {num2str(ulim(2),numSigDig), num2str(ulim(1),numSigDig)};
        {'Plant LTI'}, {'sys'}, {'fpid_optimize_lti'}
        };
    
    % Assign required workspace variables
    assignin('base', 'fpid_optimize_fpid', Gc);       % Controller
    assignin('base', 'fpid_optimize_lti', ltimodel);  % Plant model
    assignin('base', 'fpid_optimize_ref', refval);    % Reference signal
    
    % System change check
    systemChanged = 0;
    
    % Set all required model parameters if corresponding
    % blocks are present in model
    for n=1:size(PARAMS,1)
        
        % Fetch block parameters
        cParam = PARAMS(n,:);
        cParamName = cParam{1};
        cParamParams = cParam{2};
        cParamValues = cParam{3};
        
        [ig1, modelNameNoExt] = fileparts(sopt.modelname);
        cParamFullname = char(strcat(modelNameNoExt,'/',cParamName));
        
        % Check if block exists
        if ismember(cParamFullname, blks)
            % When found, set the requested block properties
            for m=1:numel(cParamParams)
                % Check if parameter not already set
                if ~strcmp(get_param(cParamFullname, cParamParams{m}), cParamValues{m})
                    set_param(cParamFullname, cParamParams{m}, cParamValues{m});
                    systemChanged = 1;
                end
            end
        end
    end
    
    % Save system, if changed
    if systemChanged
        save_system(modelName);
    end
    
    % Simulink compatibility check
    if (mver() < 7.8)
       
        simOptStruct = simset('Solver', 'ode23tb', ...
                              'MinStep', 'auto', ...
                              'MaxStep', sopt.dtmax);
                          
        [t_sim, x, y_sim, u_sim, r_sim] = sim(modelName, sopt.tmax, simOptStruct);
        y = y_sim; u = u_sim; t = t_sim; r = r_sim;
        
    else
    
        % Simulate response
        y_sim = sim(modelName, 'StartTime', '0', ...
                               'StopTime', num2str(sopt.tmax, numSigDig), ...
                               'Solver', 'ode23tb', ...
                               'MinStep', 'auto', ...
                               'MaxStep', num2str(sopt.dtmax, numSigDig));
        
        % Get response parameters
        y_sim_all = y_sim.get('fpid_optimize_y');
        y = y_sim_all(:,1);                         % System response
        u = y_sim_all(:,2);                         % Control law
        r = y_sim_all(:,3);                         % Reference

        t_sim = y_sim.get('fpid_optimize_t');       % Simulation time vector
        t = t_sim;
    
    end
    
    % Close model
    close_system(modelName);

end

