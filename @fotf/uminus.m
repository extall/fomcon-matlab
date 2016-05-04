function G=uminus(G1)
% UMINUS  Urinary minus for fractional-order transfer functions
G=fotf(G1.a,G1.na,-G1.b,G1.nb);