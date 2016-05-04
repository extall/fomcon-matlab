function G=simple(G1)
% SIMPLE Fractional-order system simplification

    [a,n]=polyuniq(G1.a,G1.na);
    G1.a=a;
    na=n;
    
    [a,n]=polyuniq(G1.b,G1.nb);
    G1.b=a;
    nb=n;
    
    % Check empty FOTF
    if isempty(nb)
        nb = 0;
    end
    
    if isempty(na)
        na = 0;
    end
    
    if isempty(G1.b)
        G1.b = 0;
    end
    
    if isempty(G1.a)
        G1.a = 0;
    end
        
    nn = min(na(end),nb(end));
    nb=nb-nn;
    na=na-nn;
    
    G = fotf(G1.a, na, G1.b, nb, G1.ioDelay);
    
end
    
function [a,an]=polyuniq(a,an)
    
    [an,ii]=sort(an,'descend');
    a=a(ii);
    ax=diff(an);
    key=1;
    
    for i=1:length(ax)
        if ax(i)==0,
            a(key)=a(key)+a(key+1);
            a(key+1)=[];
            an(key+1)=[];
        else
            key=key+1;
        end
    end
end