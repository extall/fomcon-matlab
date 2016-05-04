function [Kp,Ki,Kd] = foipdttune(K,L,T)
%FOIPDTTUNE PID tuning algorithm for FOIPDT models
%
%   Usage:  [Kp,Ki,Kd] = foipdttune(K,L,T)

    Kp = 1.111*T/(K*L^2)*1/(1+(T/L)^0.65);
    Ki = Kp / (2*L*(1+(T/L)^0.65));
    Kd = Kp * ((2*L*(1+(T/L)^0.65)) / 4);

end

