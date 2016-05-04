function G1 = mpower(G,n)
% MPOWER  Matrix power for fractional-order transfer functions

    if n == fix(n)
        
        if n>=0
            G1 = 1;
            for i=1:n
                G1=G1*G;
            end
        else
            G1=inv(G^(-n));
        end
        
        % Convert to fotf if computations
        % result in a numeric value
        if isnumeric(G1)
            G1 = fotf(G1);
        end
        
        G1 = simple(G1);
        G1.ioDelay=n*G.ioDelay;
        
    elseif length(G.a)*length(G.b)==1 && G.na==0 && G.nb==1
        
        G1=fotf(1,0,1,n);
        
    else
        error('MPOWER:NotInteger','Power must be an integer!');
    end

end

