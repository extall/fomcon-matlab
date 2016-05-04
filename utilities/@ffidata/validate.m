function varargout = validate(id, G)
%VALIDATE Validate fractional model frequency-domain identification results
%
% Usage:  [ERR_M, ERR_P] = VALIDATE(ID, G)
%         [ERR_ABS] = VALIDATE(ID, G)
%
%         where ERR_M (magnitude), ERR_P (phase) are the optional output
%                   arguments which contain the absolute error vectors
%                   G(jw) - G_id(jw) at frequency points w specified in ID.
%                   If omitted, plots validation results.
%
%               ERR_ABS is error vector containing absolute values of
%               complex response, i.e. abs(r(i))
%
%               ID is a FFIDATA structure with validation data
%
%               G is the model to validate
%
% See also: ffidata, ffid

    % Initialize output arguments
    varargout = {};

    % Number of arguments
    if nargin < 2
        error('VALIDATE:NotEnoughInputArguments','Not enough input arguments');
    end
    

    % Simulate system
    [mag, ph] = bode(G, id.w);
    [mag2, ph2] = bode(id);

    % Squeeze data
    mag = squeeze(mag);
	ph = squeeze(ph);

    mag2 = squeeze(mag2);
    ph2 = squeeze(ph2);

    % Calculate error
    err_m = mag - mag2;
    err_p = ph - ph2;


    % Check no. of output arguments
    if nargout == 0
        
        % Plot results
        h = figure();

        subplot(4,1,1);
        semilogx(id.w, mag2db(mag2), id.w, mag2db(mag), '--', 'Linewidth', 2);
        ylabel('Magnitude [dB]');
        legend('Initial data','Identified model','Location','Best');
        grid;

        subplot(4,1,2);
        semilogx(id.w, ph2, id.w, ph, '--', 'Linewidth', 2);
        ylabel('Phase [deg]');
        legend('Initial data','Identified model','Location','Best');
        grid;
        
        subplot(4,1,3);
        semilogx(id.w, err_m, 'Color', 'red');
        ylabel('Magnitude error');
        grid;
        
        subplot(4,1,4);
        semilogx(id.w, err_p, 'Color', 'red');
        xlabel('Time [s]');
        ylabel('Phase error [deg]');
        grid;

        set(h, 'NumberTitle', 'off');
        set(h, 'Name', 'Frequency-domain validation results');
        
    else
        
        varargout{1} = err_m;
        varargout{2} = err_p;
        
    end     
    
end