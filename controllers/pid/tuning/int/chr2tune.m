function [Kp,Ki,Kd] = chr2tune(K,L,T)
%CHR2TUNE Disturbance rejection with 20% overshoot CHR PID tuning formula
%  Usage: [Kp,Ki,Kd] = chr2tune(K,L,T)
%
%  Note: this function currently computes PID controller parameters,
%        P/PI computations are not implemented

    a = K*L/T;
    
    Kp = 1.2 / a;
    Ki = Kp / (2 * L);
    Kd = (0.42 * L) * Kp;

end

