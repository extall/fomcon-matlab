classdef refgen
    %REFGEN Generate a sequence of signals with the corresponding time vector
    %
    %      Usage: ST = REFGEN(...)
    %
    %      OUTPUT: ST - generated signal sequence structure with the
    %                   following parameters:
    %
    %      ST.u             % Generated signal sequence
    %      ST.t             % Generated times vector
    %      ST.signals       % Signal list
    %      ST.numSignals    % Number of signals in sequence
    %      ST.Ts            % Sampling interval
    %
    %      INPUT: series of parameters in form SIGNAL_TYPE, PARAMS, TIME which
    %             can be the following
    %      'prbs', [Tt A OFF S], T - PRBS7 sequence, Tt- single tap period [s],
    %                           A - amplitude, OFF - offset, S - 7bit starting
    %                           sequence of total duration T (alias:
    %                           'prbs7').
    %      'c', A, T - constant value of amplitude A, total duration T.
    %      'sin', [F A PHI OFF], T - sine wave of frequency F [Hz], amplitude A,
    %                                phase shift PHI and offset OFF with total
    %                                duration T.
    %      'sq', [F A PHI OFF], T - square wave signal.
    %      'tr', [F A PHI OFF], T - triangle wave signal.
    %      'swt', [F A PHI OFF], T - sawtooth wave signal.
    %      'r', [FROM TO], T - ramp signal which starts at FROM and goes to TO
    %                          over the duration T
    %      NB! Supply an empty matrix [] for any set of parameters
    %          to use default values.
    %      'Ts', T - Sampling interval T
    
    properties
        u             % Generated signal sequence
        t             % Generated times vector
        signals       % Signal list
        numSignals    % Number of signals in sequence
        Ts            % Sampling interval
    end
    
    methods
        function s = refgen(varargin)
            
            % Accepted signals
            acceptedSignals = {'prbs', 'prbs7', 'c', 'sin', 'r', 'sq', 'tr', 'swt'};
            
            % Parse input structure
            st = struct;
            st.signals=[];
            st.Ts = 0.01; % Default sample rate
            k = 1;
            c = 1;
            
            while(k<nargin)
                sig = varargin{k};
                k=k+1;
                % Sampling interval
                if (strcmpi(sig, 'ts'))
                    st.Ts = varargin{k};
                    % Particular signal
                elseif any(strcmpi(sig, acceptedSignals))
                    params = varargin{k}; k=k+1;
                    duration = varargin{k};
                    st.signals.(['s' num2str(c)]) = {lower(sig), params, duration};
                    c=c+1;
                else
                    error('REFGEN:UnknownSignalType', ...
                        'Unknown signal type specified.');
                end
                k=k+1;
            end
            
            % Number of signals
            st.N = c-1;
            
            % Order fields
            st = orderfields(st);
            
            % Time vector
            t = 0;
            
            % Signal vector
            u = [];
            
            for k=1:st.N
                
                mysig = st.signals.(['s' num2str(k)]);
                
                switch mysig{1}
                    case {'prbs7', 'prbs'}
                        
                        params = mysig{2};
                        dur = mysig{3};
                        
                        % Default parameters
                        T   = 1;
                        A   = 1;
                        OFF = 0;
                        S   = 2;
                        
                        % Check number of parameters
                        if ~isempty(params)
                            if length(params)>3, S = params(4); end
                            if length(params)>2, OFF = params(3); end
                            if length(params)>1, A = params(2); end
                            T = params(1);
                        end
                        
                        % Set the default parameters
                        st.signals.(['s' num2str(k)]) = {mysig{1}, [T A OFF S], dur};
                        
                        nN = floor(dur/T);
                        samples = A*genprbs7(nN, S);
                        numSamples = dur/st.Ts;
                        u1 = zeros(numSamples, 1);
                        sampleLength = floor(numSamples/nN);
                        
                        % Assign values
                        for n=0:nN-1
                            u1((n*sampleLength+1):(sampleLength*(n+1))) = samples(n+1);
                        end
                        
                        u = [u; OFF+u1];
                        t = [t; (t(end)+st.Ts:st.Ts:(t(end)+dur))'];
                        
                    case {'sin', 'tr', 'sq', 'swt'}
                        
                        params = mysig{2};
                        dur = mysig{3};
                        
                        % Default parameters
                        F = 1;
                        A = 1;
                        PHI = 0;
                        OFF = 0;
                        
                        % Check number of parameters
                        if ~isempty(params)
                            if length(params)>3, OFF = params(4); end
                            if length(params)>2, PHI = params(3); end
                            if length(params)>1, A = params(2); end
                            F = params(1);
                        end
                        
                        % Set the default parameters
                        st.signals.(['s' num2str(k)]) = {mysig{1}, [F A PHI OFF], dur};
                        
                        % Convert parameters
                        w  = hz2rads(F);
                        ph = deg2rad(PHI);
                        
                        % Choose appropriate function
                        switch mysig{1}
                            case 'sin'
                                g_fun = @(x) sin(x);
                            case 'tr'
                                g_fun = @(x) triangle_1(x);
                            case 'sq'
                                g_fun = @(x) square_1(x);
                            case 'swt'
                                g_fun = @(x) sawtooth_1(x);
                        end
                        
                        u1 = OFF + A*g_fun(w*(0:st.Ts:dur)'+ph);
                        u = [u; u1];
                        t = [t; (t(end)+st.Ts:st.Ts:(t(end)+dur))'];
                        
                    case 'c'
                        
                        val = mysig{2};
                        dur = mysig{3};
                        
                        % Default value
                        C = 1;
                        
                        if ~isempty(val), C = val; end
                        
                        % Set the default parameters
                        st.signals.(['s' num2str(k)]) = {mysig{1}, val, dur};
                        
                        % Generate the constant signal
                        numSamples = dur/st.Ts;
                        u1(1:numSamples,1) = C;
                        u = [u; u1];
                        t = [t; (t(end)+st.Ts:st.Ts:(t(end)+dur))'];
                        
                    case 'r'
                        
                        params = mysig{2};
                        dur = mysig{3};
                        
                        % Default value
                        FROM = 0;
                        TO = 1;
                        
                        if ~isempty(params)
                            if length(params)>1, TO = params(2); end
                            FROM = params(1);
                        end
                        
                        % Set the default parameters
                        st.signals.(['s' num2str(k)]) = {mysig{1}, [FROM TO], dur};
                        
                        % Generate ramp
                        numSamples = dur/st.Ts;
                        du = (TO-FROM)/numSamples;
                        u1 = (FROM:du:TO)';
                        u = [u; u1];
                        t = [t; (t(end)+st.Ts:st.Ts:(t(end)+dur))'];
                        
                end
                
                % Remove last sample(s) for consistency in either vector
                if length(t)>length(u), t = t(1:length(u)); end
                if length(u)>length(t), u = u(1:length(t)); end
                
            end
            
            % Assign properties
            s.u = u;
            s.t = t;
            s.signals = st.signals;
            s.numSignals = st.N;
            s.Ts = st.Ts;
            
        end
    end
    
end

% Implementation of the PRBS7 sequence
function seq = genprbs7(nper, st)

seq = zeros(nper, 1);

% Make sure we are dealing with unsigned integer
st = uint8(st);

for k=1:nper
    seq(k) = bitand(st, 1);
    nbit   = bitxor(bitand(st, 1),bitshift(bitand(st,2),-1));
    st     = bitor(bitshift(st, -1), bitshift(nbit,6));
end

st = double(st);

end

% Square wave
function x = square_1(t)
x = sign(sin(t));
end

% Sawtooth wave
function x = sawtooth_1(t)
x = 2*(t/(2*pi)-floor(1/2+t/(2*pi)));
end

% Triangle wave
function x = triangle_1(t)
x = 2*abs(sawtooth_1(t))-1;
end


