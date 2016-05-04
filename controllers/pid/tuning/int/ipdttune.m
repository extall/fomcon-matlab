function [Kp,Ki,Kd] = ipdttune(K,L)
%IPDTTUNE IPDT model ISE based PID tuning
%  Usage: [Kp,Ki,Kd] = ipdttune(K,L)

    Kp = 1.37 / (K*L);
    Ki = Kp / (1.49 * L);
    Kd = Kp * (0.59 * L);

end

