function G = fogpm(K,L,T,type,alpha)
% FOGPM Get FO process model specified by gain, delay,time constant,order.
% Usage: FOGPM(K,L,T,type,alpha) Returns a plant model, specified by type:
%                        'ffopdt'    G = K / (1+T*s^alpha) * exp(-L)
%                        'fipdt'     G = K / s^alpha * exp(-L)
%                        'ffoipdt'   G = K / s^alpha * (1+T*s) * exp(-L)

    % Check input arguments
    if nargin < 5
        alpha = 1;
    end
    
    s = fotf('s')^alpha;

    switch lower(type)
        
        case 'ffopdt'

            G = K / (1+T*s);
            G.ioDelay = L;
            
        case 'fipdt'
            
            G = K / s;
            G.ioDelay = L;
            
        case 'ffoipdt'
            
            G = K / (fotf('s')*(1+T*s));
            G.ioDelay = L;
            
    end
    
end

