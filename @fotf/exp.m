function G = exp(H)
%EXP Delay for fractional-order transfer functions
    
    % Use MATLAB built-in function to evaluate the delay
    s=tf('s');
    [a,na,b,nb]=fotfparam(H);
    to_eval = fpoly2str(b,nb,'s');
    
    % Check for +/-"1"s, i.e. for +s and -s cases
    if ~strcmpi(to_eval, '-s') && ~strcmpi(to_eval, '+s')
       to_eval = strrep(to_eval,'s','*s');
    end
    
    % Check zero input value
    if ~isempty(to_eval)
        % Proceed with evaluation
        to_eval = ['exp(' to_eval ')'];
        Z = eval(to_eval);
    else
        Z.outputDelay = 0;
    end

    G = fotf(1,0,1,0,Z.outputDelay);
    
end

