function varargout = frac_struct_rc_cauer2(R,C,~,K,params)
%FRAC_STRUCT_RCL_CAUER2 Impedance of Cauer 2 type RC network

    % Standard output scheme: return type of network in form [R C L]
    if nargin == 0
        varargout{1} = [1 1 0];
        return;
    end
    
    % Check that the length of input vectors
    if length(R) < 1 || length(C) < 1
        error('FracStructRcCauer2:InsufficientNumberOfComponents', ...
        'Insufficient number of components provided');
    end
    
    if length(R) ~= length(C)
        error('FracStructRcCauer2:ComponentVectorsNotOfSameLength', ...
        'Component vectors must be of the same length');
    end
    
    % Initial function
    Z = zpk((1/C(end))/zpk('s') + 1/(1/R(end)));
    
    for n=length(R)-1:-1:1
        Z = 1/C(n)/zpk('s') + 1/(1/R(n) + 1/Z);
    end
    
    % Return the impedance Z(s)
    varargout{1} = K * Z;
    
end

