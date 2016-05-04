function [Kp,Ki,lambda] = fmigotune(K,L,T)
%FMIGOTUNE FMIGO PI^lambda tuning formula
%   Usage: [Kp,Ki,lambda] = fmigotune(K,L,T)

    % Check input arguments
    if nargin < 3
        error('FMIGOTUNE:NotEnoughInputArguments', ...
              'Not enough input arguments');
    end

    % Relative dead time
    tau = L / (L+T);

    % Determine lambda
    if (tau >= 0.6)
        lambda = 1.1;
    elseif (0.4 <= tau && tau < 0.6)
        lambda = 1.0;
    elseif (0.1 <= tau && tau <0.4)
        lambda = 0.9;
    elseif (tau < 0.1)
        lambda = 0.7;
    end

    % Determine controller gains
    Kp = 1/K * (0.2978/(tau+0.00307));
    Ki = Kp / (T*(0.8578/(tau^2-3.402*tau+2.405)));

end

