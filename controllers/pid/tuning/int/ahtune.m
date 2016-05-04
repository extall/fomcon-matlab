function [Kp,Ki,Kd] = ahtune(K,L,T)
%AHTUNE Astrom-Hagglund 'AMIGO' PID tuning formula
%    Usage: [Kp,Ki,Kd] = ahtune(K,L,T)
%
%    Note: This function currently computes PID controller parameters

    Kp = 1/K*(0.2+0.45*T/L);
    Ki = Kp/((0.4*L+0.8*T)/(L+0.1*T)*L);
    Kd = (0.5*L*T/(0.3*L+T))*Kp;

end

