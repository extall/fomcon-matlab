function [Kp, Ki, Kd] = cctune(K,L,T)
%CCTUNE Cohen-Coon PID tuning formula for PID controllers
%  Usage: [Kp,Ki,Kd] = cctune(K,L,T)
%
%  Note: this function currently computes PID controller parameters,
%        P/PI/PD computations are not implemented

    a = K * L / T;
    tau = L / (L + T);

    Kp = (1.35/a)*(1+(0.18*tau)/(1-tau));
    Ki = Kp/((2.5-2*tau)/(1-0.39*tau)*L);
    Kd = ((0.37-0.37*tau)/(1-0.81*tau)*L)*Kp;

end

