function [Kp, Ki, Kd] = chr1tune(K, L, T)
%CHR1TUNE Set-point regulation with 20% overshoot CHR PID tuning formula
%  Usage: [Kp,Ki,Kd] = chr1tune(K,L,T)
%
%  Note: this function currently computes PID controller parameters,
%        P/PI computations are not implemented

    a = K*L/T;
    
    Kp = 0.95 / a;
    Ki = Kp/(1.4 * T);
    Kd = (0.47 * L)*Kp;

end

