function G=minus(G1,G2)
% MINUS  Subtraction for fractional-order input-output models
    if strcmp(class(G1), 'double')
        newFotf = fotf(1,0,G1,0);
        G1      = newFotf;
    end
    
    if strcmp(class(G2), 'double')
        newFotf = fotf(1,0,G2,0);
        G2      = newFotf;
    end

    G=G1+(-G2);
end