function G = new_fod(r, N, wb, wh, b, d)
%NEW_FOD Creates a new FO transfer function approximation in ZPK
%
%   NEW_FOD(r, N, wb, wh, b, d) creates a ZPK object with a refined
%          Oustaloup filter approximation of a fractional-order operator
%          s^r of order N and valid within frequency range (wb, wh).
%
%          s^r = (d*wh/b)^r * (ds^2 + b*wh*s)/(d*(1-r)*s^2+b*wh*s+d*r)*Gp,
%          where       
%                      N
%                    -----
%               Gp =  | |  (s+w'k)/(s+wk),
%                     | |
%                    k= -N
%
%               wk  = (b*wh/d)^((r+2k)/(2N+1)),
%               w'k = (d*wb/b)^((r-2k)/(2N+1)).
%
%          Should parameters b and d be omitted, the default values will be
%          used: b=10, d=9.
	
	% Number of input arguments control
	if nargin < 4
		error('NEW_FOD:NotEnoughInputArguments', 'Not enough input arguments.');
	end
	
    % Check input arguments, if b and d are
    % omitted, use default (recommended) values
    if nargin==4
        b=10;
        d=9;
    end
    
    % If r=0, return empty zpk
    if r == 0
        
        G = zpk;
        
    else

        mu = wh/wb;
        k=-N:N;
        w_kp=(mu).^((k+N+0.5-0.5*r)/(2*N+1))*wb;
        w_k=(mu).^((k+N+0.5+0.5*r)/(2*N+1))*wb;
        K=(d*wh/b)^r;
        G=zpk(-w_kp',-w_k',K)*tf([d,b*wh,0], [d*(1-r),b*wh,d*r]);
        
    end

end
