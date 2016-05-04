function [Kp,Ki,Kd] = zntune(K,L,T,ctrlType)
%ZNTUNE Ziegler-Nichols tuning algorithm based on time response
%
%       Usage: [Kp,Ki,Kd] = ZNTUNE(K,L,T,ctrlType)
%       
%       where   [Kp, Ki, Kd] - controller parameters,
%               K, L, T      - FOPDT model parameters
%               ctrlType     - 1, 2 or 3 for P, PI or PID controller
    
    % Initial values
    Kp = 0;
    Ki = 0;
    Kd = 0;
    
    % Get parameters
    a = K * L / T;
    
    % Return parameters based on the controller type
    switch ctrlType
        
        case 1
            % P controller
            Kp = 1/a;
            
        case 2
            % PI controller
            Kp = 0.9/a;
            Ki = Kp / (3 * L);
            
        case 3
            % PID controller
            Kp = 1.2/a;
            Ki = Kp / (2 * L);
            Kd = (L / 2) * Kp;
            
        otherwise
            % Wrong switch
            error('Invalid controller type!');    
    end

end

