function Gc = tid(Kt, n, Ki, Kd)
%TID Tilt-Integral-Derivative controller
%
%    This function will return a TID controller in FOTF form.
%
%    Usage: GC = TID(KT, N, KI, KD)
%
%    A controller is obtained with the following structure:
%    GC = KT / s^(1/N) + KI/s + KD*s

    % Check number of parameters
    if nargin < 4
        error('TID:NotEnoughInputArguments', ...
              'Not enough input arguments.');
    end
    
    Gc  =   fotf(1,1/n,Kt,0) + ...
            fotf(1,1,Ki,0) + ...
            fotf(1,0,Kd,1);

end

