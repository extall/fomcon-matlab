function z = ipdtfun(x, y, t)
%IPDTFUN Helper function for IPDT fitting
%    USAGE: z=ipdtfun(x, y, t)

    % Get parameters  
    K = x(1);
    L = x(2);
    
    % Build model
    G = gpm(K,L,0,'ipdt');
    
    % Get time response
    y1 = step(G, t);
    
    % Get error
    z = y-y1;

end

