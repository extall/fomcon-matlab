function varargout = frac_struct_rc_lowpass_c1(R,C,~,K,params)
%FRAC_STRUCT_RC_LOWPASS Impedance of a fractional low-pass filter

    % -----STRUCTURE USES -------------------------------------------------
    STRUCTURE1 = 'frac_struct_rc_cauer1';
    % ---------------------------------------------------------------------

    % Standard output scheme: return type of network in form [R C L]
    if nargin == 0
        varargout{1} = [1 1 0];
        return;
    end
    
    % Check that the length of input vectors
    if length(R) < 2 || length(C) < 2
        error('FracStructRcCauer1:InsufficientNumberOfComponents', ...
        'Insufficient number of components provided');
    end
    
    % Get first resistor value
    R1 = R{1};
    R  = R{2};
    C  = C{2};
    K1 = K{1};
    K2 = K{2};
    
    % Get Cauer 1 network impedance
    str_fun = str2func(STRUCTURE1);
    Z = 1 / (1 + R1 * 1/str_fun(R,C,[],K2));
    varargout{1} = K1 * Z;

end

