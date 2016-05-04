function G = mrdivide(G1,G2)
% MRDIVIDE Right division for fractional-order dynamic systems.
%
% Note: if delays in the systems are present, dividing the two systems may
% result in positive delays in which case the overall delay of the system
% will be changed to zero with a warning.

    G1=fotf(G1);
    G2=fotf(G2);
    G=G1*inv(G2);
    G.ioDelay=G1.ioDelay-G2.ioDelay;
    
    if G.ioDelay<0
        G.ioDelay = 0;
        warning('MRDIVIDE:PositiveDelay','Resulting block has positive delay: changing to zero');
    end
end

