function varargout = frac_struct_rc_foster1(R,C,~,K,params)
%FRAC_STRUCT_RCL_FOSTER1 Impedance of Foster 1 type RC network

    % Standard output scheme: return type of network in form [R C L]
    if nargin == 0
        varargout{1} = [1 1 0];
        return;
    end
    
    % Check that the length of input vectors
    if length(R) < 2 || length(C) < 2
        error('FracStructRcFoster1:InsufficientNumberOfComponents', ...
        'Insufficient number of components provided');
    end
    
    if length(R) ~= length(C)
        error('FracStructRcFoster1:ComponentVectorsNotOfSameLength', ...
        'Component vectors must be of the same length');
    end
    
    % Form the impedance Z(s)
    Z = zpk(R(1)) + 1/C(1)/zpk('s');
    
    for n=2:length(C)
        K_n = 1/C(n);
        sigma_n = K_n/R(n);
        Z = Z + K_n/(zpk('s')+sigma_n);
    end
    
    % Return the impedance Z(s)
    varargout{1} = K * Z;
    
end

