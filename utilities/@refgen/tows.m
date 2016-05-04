function st = tows(r,op)
%TOWS Convert REFGEN object to structure as needed by the From Workspace
%     block in Simulink.
%
% Usage: TOWS(R,OP) where R is the REFGEN object, and
%                         OP are additional options.
% Currently supported options:
%            OP.time = 't'   Populate the time vector
%                    = 'e'   Leave time vector empty (default)

% Set default values for options
if nargin<2
	gt = 0;			% Do not populate the time vector by default
else
	% Go through options
    if cfieldexists(op, 'time')
        switch(op.time)
            case 't'
                gt = 1;
            case 'e'
                gt = 0;
        end
    end
end



st = struct;

% Populate the time vector
if (gt)
    st.time = r.t;
else
    st.time = [];
end

st.signals.values = r.u;
st.signals.dimensions = 1;
st.ts = r.Ts;

end

