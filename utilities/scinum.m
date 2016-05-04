function y = scinum(x)
%SCINUM Convert a floating point number to scientific format A * 10^B

    % Check input vector
    [rows, cols] = size(x);
    if rows > 1 && cols > 1
        error('SCINUM:NotAVector', 'Input argument is not a vector.');
    end
    
    % Convert to column vector form, if needed
    if rows == 1 && cols > 1, x = x'; end
    
    exponent = floor(log10(abs(x)));
    y = [x./(10.^exponent) exponent];

end

