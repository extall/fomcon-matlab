function G=plus(G1,G2)
% PLUS  Adds two fractional-order input/output models together.
%
% Note, that you cannot add two systems G1 and G2 if they exhibit
% different I/O delays.

    G1 = fotf(G1);
    G2 = fotf(G2);
    
    na=[];
    nb=[];
    
    if G1.ioDelay == G2.ioDelay
        a=kron(G1.a,G2.a);
        b=[kron(G1.a,G2.b),kron(G1.b,G2.a)];

        for i=1:length(G1.a)
            na=[na G1.na(i)+G2.na];
            nb=[nb, G1.na(i)+G2.nb];
        end

        for i=1:length(G1.b)
            nb=[nb G1.nb(i)+G2.na];
        end

        G=simple(fotf(a,na,b,nb,G1.ioDelay));
    else
       error('PLUS:DifferentDelays', 'Cannot handle different delays'); 
    end
    
end