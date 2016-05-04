function G=mtimes(G1,G2)
% Multiplies two fractional-order input/output models together.
    
    G1 = fotf(G1);
    G2 = fotf(G2);
    
    na=[];
    nb=[];

    a=kron(G1.a,G2.a);
    b=kron(G1.b,G2.b);
    
    for i=1:length(G1.na)
        na=[na,G1.na(i)+G2.na];
    end

    for i=1:length(G1.nb)
        nb=[nb,G1.nb(i)+G2.nb];
    end

    G=simple(fotf(a,na,b,nb,G1.ioDelay+G2.ioDelay));
    
end