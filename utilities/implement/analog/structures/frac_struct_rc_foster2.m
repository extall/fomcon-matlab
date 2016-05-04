function varargout = frac_struct_rc_foster2(R,C,~,K,params)
%FRAC_STRUCT_RCL_FOSTER2 Impedance of Foster 2 type RC network

    % Standard output scheme: return type of network in form [R C L]
    if nargin == 0
        varargout{1} = [1 1 0];
        return;
    end
    
    % Check that the length of input vectors
    if length(R) < 2 || length(C) < 2
        error('FracStructRcFoster2:InsufficientNumberOfComponents', ...
        'Insufficient number of components provided');
    end
    
    if length(R) ~= length(C)
        error('FracStructRcFoster2:ComponentVectorsNotOfSameLength', ...
        'Component vectors must be of the same length');
    end
    
    % Form the admittance Y(s)
    Y = zpk(1/R(1)) + C(1)*zpk('s');
    
    for n=2:length(R)
        K_n = 1/R(n);
        sigma_n = K_n/C(n);
        Y = Y + K_n*zpk('s')/(zpk('s')+sigma_n);
    end
    
    % Return the impedance Z(s) = 1/Y(s)
    varargout{1} = K * (1 / Y);
    
end

