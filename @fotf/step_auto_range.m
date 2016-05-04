function t = step_auto_range(G)
%STEP_AUTO_RANGE Automatically locate appropriate step response length

% No. of points for auto-ranging (if dcgain exists) and vector generation
AUTORANGE_NUM_PTS_AUTO = 10;
AUTORANGE_NUM_PTS_GEN  = 10000;

% Detection coefficient "final value is at least 90% of dc gain"
AUTORANGE_DET_COEF = 0.9;

% Add this value to the exponent in case of oscillating processes
AUTORANGE_OSC_ADD  = 1;

% Default end point in case dcgain gives "0" or "Inf"
AUTORANGE_DEFAULT_END = 100;

% Test DC gain
myGain = dcgain(G);
if isinf(myGain) || abs(myGain)<eps
    t = linspace(0,AUTORANGE_DEFAULT_END,AUTORANGE_NUM_PTS_GEN);
else
    % Final time range locator
    t_exp = -4;
    t_exp_max = 9;
    while t_exp<t_exp_max
        t = linspace(0, 10^t_exp, AUTORANGE_NUM_PTS_AUTO);
        u = ones(size(t));
        y = lsim(G,u,t);
        if abs(y(end)/abs(myGain))>=AUTORANGE_DET_COEF
            % Check if this is an oscillating process, and add a few
            % exponent values to ensure that it is covered in full
            if max(abs(y)) > abs(y(end))
                t_exp = t_exp + AUTORANGE_OSC_ADD;
            end
            
            % All done.
            break;
        end
        
        t_exp = t_exp + 1;
    end
    
    % Use the obtained value to generate the vector
    t = linspace(0,10^t_exp, AUTORANGE_NUM_PTS_GEN);
    
end

end

