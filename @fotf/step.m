function varargout=step(G,t)
% STEP  Step response of fractional-order dynamic systems.
%
% Usage: [Y, T] = step(G)
%
% where Y is the system response, T is the time vector, and G is the
% fractional-order SISO system.

varargout = {};

if nargin < 2
    t = step_auto_range(G);
end

u=ones(size(t));
y=lsim(G,u,t);

if size(y,2) > size(y,1)
    y = y';
end

switch (nargout)
    case 0
        plot(t,y);
        title('Step response');
        xlabel('Time [s]');
        ylabel('Amplitude');
        
        % Plot final value if present
        % Test DC gain
        myGain = dcgain(G);
        if isinf(myGain) || abs(myGain)<eps
            % Do nothing
            
        else
            % Plot final value
            hold on; plot([t(1) t(end)], [myGain myGain], ':k');
            
        end
    case 1
        varargout{1} = y;
    case 2
        varargout{1} = y;
        varargout{2} = t;
end

end
