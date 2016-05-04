function G = gpm(K,L,T,type)
% GPM Get process model specified by gain, delay and time constant.
% Usage: GPM(K,L,T,type) Returns a plant model, specified by type:
%                        'fopdt'    G = K / (1+T*s) * exp(-L)
%                        'ipdt'     G = K / s * exp(-L)
%                        'foipdt'   G = K / s * (1+T*s) * exp(-L)

    s = tf('s');

    switch lower(type)
        
        case 'fopdt'

            G = K / (1+T*s);
            G.ioDelay = L;
            
        case 'ipdt'
            
            G = K / s;
            G.ioDelay = L;
            
        case 'foipdt'
            
            G = K / (s*(1+T*s));
            G.ioDelay = L;
            
    end
    
end

