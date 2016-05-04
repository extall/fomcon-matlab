function varargout=impulse(G,t)
% IMPULSE Impulse response of fractional-order dynamic systems.
%
%   Usage: IMPULSE(G, T)
%          IMPULSE(G)
%          Y = IMPULSE(G, T)
%   where
%          G is the LTI dynamic system in FOTF form,
%          T is the regularly spaced time vector,
%          Y is the obtained response; if omitted, the impulse response of
%          the system is plotted.

if nargin < 2
    % Use the same auto-ranging function as step
    t=step_auto_range(G);
end

varargout = {};

% Compute du/dt, replace first element, and pad to the right
y = diff(step(G,t))/(t(2)-t(1)); y(1)=y(2); y(end+1)= y(end);

if nargout==0
    plot(t,y);
    title('Impulse response');
    xlabel('Time [s]');
    ylabel('Amplitude');
    
    % Plot final value if present
    % Test DC gain
    myGain = dcgain(G);
    if isinf(myGain) || abs(myGain)<eps
        % Do nothing
        
    else
        % Plot final value
        hold on; plot([t(1) t(end)], [0 0], ':k');
        
    end
    
else
    varargout{1} = y;
end

end