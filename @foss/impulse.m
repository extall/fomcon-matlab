function varargout = impulse(S, t)
% IMPULSE  Impulse response of fractional-order dynamic systems.
%
% The impulse response of multi-input systems is the collection of
% impulse responses for each input channel.
%
% See also: fotf/impulse


varargout = {};

noOuts = size(S.fotf_array,1);
noIns  = size(S.fotf_array,2);

if nargin < 2
    % Go through all of the FO TFs in the FOSS array, and find the correct
    % time range (maximum among all possible ranges)
    tms_max = 0;
    t= [];
    for k=1:noIns
        for l=1:noOuts
            tms_my = step_auto_range(getfotf(S,l,k));
            if max(tms_my) > tms_max
                tms_max = tms_my;
                t = tms_my;
            end
        end
    end
end

if nargout < 1
    
    % Only draw the plot
    
    % Plot counter
    p = 1;
    for k=1:noOuts
        for l=1:noIns
            subplot(noOuts, noIns, p);
            impulse(getfotf(S,k,l),t);
            if k==1
                title(['From: In(' num2str(l) ')'], 'FontSize', 8);
            else
                title('');
            end
            if l==1
               ylabel(['To: Out(' num2str(k) ')'], 'FontSize', 8);
            else
               ylabel('');
            end
            xlabel('Time [s]', 'FontSize', 8);
            p=p+1;
        end
    end

elseif nargout > 1
    
    % Output array
    y = zeros(numel(t),noOuts,noIns);
    
    % Simulate all systems
    for k=1:noIns
        for l=1:noOuts
            y(:,l,k) = impulse(getfotf(S,l,k),t);
        end
    end
    
    switch (nargout)
        case 1
            varargout{1} = y;
        case 2
            varargout{1} = y;
            varargout{2} = t;
    end
    
end


end

