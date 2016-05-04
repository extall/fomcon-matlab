function G=feedback(F,H)
% FEEDBACK Feedback connection of two input/output fractional-order transfer functions. 
%
%    M = feedback(G,H) computes a closed-loop model C for the feedback loop:
%
%          u --->O---->[ G ]---------> y
%                ^               |
%                |               |
%                ------[ H ]<-----
%
%    Note, that fractional-order systems in the feedback loop must have
%    equal delays or lack thereof.
    
    F=fotf(F);
    H=fotf(H);
    na=[];
    nb=[];
    
    if F.ioDelay == H.ioDelay
    
        b=kron(F.b,H.a);
        a=[kron(F.b,H.b), kron(F.a,H.a)];

        for i=1:length(F.b)
            nb=[nb, F.nb(i)+H.na];
            na=[na, F.nb(i)+H.nb];
        end

        for i=1:length(F.a)
            na=[na F.na(i)+H.na];
        end
    
        G = simple(fotf(a,na,b,nb,F.ioDelay));
        
    else
        error('FEEDBACK:DifferentDelays','Cannot handle different delays'); 
    end
end