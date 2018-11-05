classdef ffidata
    %FFIDATA Frequency-domain identification parameters for FOTF systems
    %
    % Usage:  P = FFIDATA(MAG, PHASE, W)
    %         or
    %         P = FFIDATA(R, W)
    %         or
    %         P = FFIDATA(FRD)
    %
    % where
    %         MAG - observed frequency response magnitude in dB at W
    %         PHASE - observed frequency response phase in deg at W
    %         R - (alternative) - complex response of the system at W
    %         W - frequency vector where the response is known [rad/s] 
    %
    % Note:   If parameters are not equally sized (i.e. if parameters are
    %         obtained by means of bode function), the corresponding
    %         responses will be squeezed.
    %
    % See also: freqresp, frd, bode, levy, vinagre, hartley, ffidata/bode,
    %           ffidata/validate
    
    properties
        mag       % Observed frequency response magnitude in dB at W
        phase     % Observed frequency response phase in deg at W
        w         % Frequency vector where the response is known [rad/s] 
		tstmp     % Timestamp showing date/time of dataset creation
        focus     % Focus on complex response, mag only, or phase only
    end
    
    methods
        
        % Initializer
        function p = ffidata(varargin)
           
            % Check number of arguments
            if size(varargin,2) == 1
                frdobj = varargin{1};
                w = get(frdobj, 'Frequency');
                [mag, phase] = bode(frdobj,w);
                mag = mag2db(mag);
            elseif size(varargin,2) == 2
                w = varargin{2};
                [mag, phase] = bode(frd(varargin{1},w),w);
                mag = mag2db(mag);
            elseif size(varargin,2) == 3
                mag = varargin{1};
                phase = varargin{2};
                w = varargin{3};
            else
                 error('FFIDATA:WrongNumberOfInputArguments', 'Wrong number of input arguments.');
            end
            
            % Check magnitude
            mag = squeeze(mag);
            
            if size(mag, 1) < size(mag, 2)
                mag = mag';
            end
            
            % Check phase
            phase = squeeze(phase);
            
            if size(phase, 1) < size(phase, 2)
                phase = phase';
            end
            
            % Check w
            if size(w, 1) < size(w, 2)
                w = w';
            end
            
            % Check size match
            if ~(max(size(mag)) == max(size(phase)) && ...
                    max(size(phase)) == max(size(w)))
                error('FFIDATA:VectorsNotEquallySized', 'Data vectors are not equally sized.');
            end
            
            % Set parameters
            p.mag = mag;
            p.phase = phase;
            p.w = w;
            
			% Timestamp
            p.tstmp = datestr(now);
            
            % Focus
            p.focus = 'complex';
        end
        
    end
    
end

