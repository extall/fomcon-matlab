function z = fopdtfun(x, y, t)
%FOPDTFUN Helper function for FOPDT fitting
%   USAGE:z=fopdtfun(x, y, t)

    % Get parameters  
    K = x(1);
    L = x(2);
    T = x(3);
    
    % Build model
    G = gpm(K,L,T,'fopdt');
    
    % Get time response
    y1 = step(G, t);
    
    % Get error
    z = y-y1;

end

