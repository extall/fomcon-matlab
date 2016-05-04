function G = exp(H)
%EXP Delay for fractional-order transfer functions
    
    % Use MATLAB built-in function to evaluate the delay
    s=tf('s');
    [a,na,b,nb]=fotfparam(H);
    to_eval = poly2str(b,nb,'s');
    to_eval = ['exp(' strrep(to_eval,'s','*s') ')'];
    Z = eval(to_eval);
    G = fotf(1,0,1,0,Z.outputDelay);

end

