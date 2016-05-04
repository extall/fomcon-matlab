function G=inv(G1)
% INV Fractional-order transfer function inverse
    G=fotf(G1.b,G1.nb,G1.a,G1.na,-G1.ioDelay);
end