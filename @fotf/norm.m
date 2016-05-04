function n = norm(G, eps0)
%NORM Computes norms of given FOTF system
%
%   Usage: N = NORM(G, EPS0)
%
%   where EPS0 is either an error tolerance for computation of the H_2
%   norm, or Inf for computation of the H_inf norm.

	config = fomcon('config');
    epsi   = config.Core.General.Internal_computation_accuracy;

    dx = 1;
    f0 = 0;
    
    % Check input arguments
    if nargin == 1
        % H_2 norm computation
        eps0 = epsi;
    end
    
    % H_inf norm computation
    if nargin == 2 && ~isfinite(eps0)
       f = @(x)(-abs(freqresp(G, 1i*x)));
       x = -fminsearch(f,0);
       n = abs(freqresp(G, 1i*x));
    else
       % H_2 norm computation
       f = @(x)freqresp(G, x).*freqresp(G, -x);
       while (1)
           n = quadgk(f, -dx*1i, dx*1i) / (2*pi*1i);
           if abs(n-f0) < epsi
               break;
           else
               f0=n;
               dx = dx * 1.2;
           end
       end
    end
end

