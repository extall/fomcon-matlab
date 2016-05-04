function G1 = normalize(G)
%NORMALIZE Fractional-order transfer function coefficient factoring

    [a,na,b,nb,L] = fotfparam(G);
    factor = a(end);
    a = a ./ factor;
    b = b ./ factor;
    G1 = fotf(a,na,b,nb,L);

end

