function varargout = frac_struct_rc_cauer1(R,C,~,K,params)
%FRAC_STRUCT_RCL_CAUER1 Impedance of Cauer 1 type RC network

    % Standard output scheme: return type of network in form [R C L]
    if nargin == 0
        varargout{1} = [1 1 0];
        return;
    end
    
    % Check that the length of input vectors
    if length(R) < 1 || length(C) < 1
        error('FracStructRcCauer1:InsufficientNumberOfComponents', ...
        'Insufficient number of components provided');
    end
    
    % Initial function
    if length(R) == length(C)
        Z = R(end) + 1/C(end)/zpk('s');
    elseif length(R) > length(C)
        Z = R(end);   
    end
    
    for n=length(R)-1:-1:1
        Z = R(n) + 1/(C(n)*zpk('s')+ 1/Z);
    end
    
    % Return the impedance Z(s)
    varargout{1} = K * Z;
    
end

