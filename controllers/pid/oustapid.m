function Gc = oustapid(Kp, Ki, lam, Kd, mu, wb, wh, N, type, red)
%OUSTAPID Obtain Oustaloup approximation of the fractional-order PID controller
%
%         This function generates an Oustaloup filter approximation of a
%         fractional PID controller with special consideration for integrator
%         and differentiator realization.
%
%  Usage: GC = OUSTAPID(KP, KI, LAM, KD, MU, WB, WH, N, TYPE, RED)
%
%  where  KP, KI, LAM, KD, MU - fractional PID parameters in
%                       
%                        Gc = Kp + Ki/s^lam + Kd*s^mu
%
%         and WB, WH, N and TYPE (optional) are Oustaloup filter
%         approximation parameters (defaults will be used if these
%         parameters are omitted).
%
%         It is also possible to invoke model reduction via BALRED()
%         function of the Control System toolbox to obtain a lower order of
%         the controller.
%
%  See also: fpid, oustapp

    if nargin == 5
        Gc = Kp + ...
		     Ki * oustapp(fotf(1,ceil(lam),1,ceil(lam)-lam)) + ...
             Kd * oustapp(fotf(1,0,1,mu));
    elseif nargin >= 9
		Gc = Kp + ...
		     Ki * oustapp(fotf(1,ceil(lam),1,ceil(lam)-lam), wb, wh, N, type) + ...
             Kd * oustapp(fotf(1,0,1,mu), wb, wh, N, type);
    else
        error('OUSTAPID:WrongNumberOfArguments', 'Wrong number of arguments.');
    end
    
    % Check if reduction via BALRED() is requested
    if nargin == 10
        if round(red) < 1
            error('OUSTAPID:CannotReducePastFirstOrder', ...
                'Model reduction failed due to too low order requested');
        end
        
        Gc = balred(Gc, round(red));
    end
    
end

