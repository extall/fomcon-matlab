function [y, str] = engnum(x)
%ENGNUM Convert a floating point number to engineering format

    % Engineering exponent values and prefixes
    ENG      = (-24:3:24)';
    ENG_PREF = {'y', 'z', 'a', 'f', 'p', 'n', 'u', 'm', '', ...
                'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'};

    % Check input vector
    [rows, cols] = size(x);
    if rows > 1 && cols > 1
        error('EngNum:NotAVector', 'Input argument is not a vector.');
    end
    
    % Convert to column vector form, if needed
    if rows == 1 && cols > 1, x = x'; end
    
    % Fetch numbers in scientific format
    x = scinum(x);
    
    % Find closest exponents
    x_exp   = x(:,2);
    str   = {};                     % Strings with values
    
    for n=1:numel(x_exp)
        x_exp_d  = x_exp(n)-ENG;
        [ig1, ind] = min(abs(x_exp_d-1));
        x_exp(n) = ENG(ind);
        x(n,1)   = x(n,1)*10^x_exp_d(ind);
    end
    
    % Replace in original matrix
    x(:,2) = x_exp;
    
    % Form string
    for n=1:numel(x_exp)   
        str{n,1} = [num2str(x(n,1)),' ',ENG_PREF{x(n,2)==ENG}];
    end
    
    y = x;

end

