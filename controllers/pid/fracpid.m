function [pidCtrl, ctrlType, ltxString] = fracpid(varargin)
%FRACPID Create a new fractional-order PID controller in parallel form.
%
%        A parallel form of the fractional-order PID controller is assumed,
%        i.e. Gc(s) = Kp + Ki/s^lambda + Kd * s^mu.
%
%   Usage: [PIDCTRL, CTRLTYPE] = FRACPID(KP, KI, LAMBDA, KD, MU)
%           CTRLTYPE will contain a string with controller type.
%   Note: no parameter may be omitted! Pass Ki and/or Kd = 0 with any lambda
%         and mu to obtain various controller types (P, PI, PD, PID).

    % Check input argument
    myvar = varargin{1};
    if ischar(myvar)
        % Load the parameters from a given file
        myParams = load(myvar);
        Kp       = str2num(myParams.FPID_Optimizer_GUI_config.FPIDParams.Kp);
        Ki       = str2num(myParams.FPID_Optimizer_GUI_config.FPIDParams.Ki);
        ilambda  = str2num(myParams.FPID_Optimizer_GUI_config.FPIDParams.Lam);
        Kd       = str2num(myParams.FPID_Optimizer_GUI_config.FPIDParams.Kd);
        dmu      = str2num(myParams.FPID_Optimizer_GUI_config.FPIDParams.Mu);
    else
        % Otherwise check number of parameters
        if nargin < 5
            error('FRACPID:NotEnoughInputArguments', 'Not enough input arguments.');
        else
            [Kp, Ki, ilambda, Kd, dmu] = deal(varargin{:});
        end
        
    end
    
    % Determine controller type and create controller
    if fleq(Ki,0) && fleq(Kd,0)
        % P controller
        ctrlType = 'P';
        ltxString = 'P';
        pidCtrl  = fotf(1,0,Kp,0);
    elseif fleq(Ki,0) && ~fleq(Kd,0)
        % PD-micro controller
        ctrlType = 'PD';
        ltxString = 'PD^{\mu}';
        pidCtrl  = fotf(1,0,Kp,0) + ...
                   fotf(1,0,Kd,dmu);
    elseif ~fleq(Ki,0) && fleq(Kd,0)
        % PI-lambda controller
        ctrlType = 'PI';
        ltxString = 'PI^{\lambda}';
        pidCtrl  = fotf(1,0,Kp,0) + ...
                   fotf(1,ilambda,Ki,0);
    elseif ~fleq(Ki,0) && ~fleq(Kd,0)
        % Full fractional controller
        ctrlType = 'PID';
        ltxString = 'PI^{\lambda}D^{\mu}';
        pidCtrl  = fotf(1,0,Kp,0) + ...
                   fotf(1,ilambda,Ki,0) + ...
                   fotf(1,0,Kd,dmu);
    end
    
end
