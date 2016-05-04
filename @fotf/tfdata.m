function [num, den, q] = tfdata(G)
%TFDATA Fractional-order transfer function data
%
% Returns numerator and denominator of a fractional transfer function in a
% rational transfer function form with commensurate order q
%
% Usage: [NUM, DEN, Q] = TFDATA(G)
%
% Note, that vectors of 'double' type are returned.
    
    % Get commensurate order and max order
    q = comm_order(G);
    
    % Get elements positions
    a1=fix_s(G.na/q);
    b1=fix_s(G.nb/q);
    
    % Numerator
    num(b1+1) = G.b;
    num = num(end:-1:1);
    
    % Denominator
    den(a1+1) = G.a;
    den = den(end:-1:1);

end

