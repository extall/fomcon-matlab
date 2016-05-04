function z = foipdtfun(x, y, t)
%FOIPDTFUN Helper function for IPDT fitting
%
%  Usage: Z = foipdtfun(X, Y, T)

    % Get parameters  
    K = x(1);
    L = x(2);
    T = x(3);
    
    % Build model
    G = gpm(K,L,T,'foipdt');
    
    % Get time response
    y1 = step(G, t);
    
    % Get error
    z = y-y1;

end

