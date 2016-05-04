function G=oustafod(r,N,wb,wh)
% OUSTAFOD Oustaloup's approximation of a s^r operator
%
%  OUSTAFOD(r,N,wb,wh) computes the Oustaloup filter approximation of a
%           fractional-order operator s^r of the order N and valid in the
%           frequency range of (wb, wh). The function returns a ZPK object
%           containing the continuous-time filter. The following equation
%           is used to construct the filter:
%                       N
%                     -----
%           s^r = K *  | | (s+w'k)/(s+wk),
%                      | |
%                     k= -N
%
%           where w'k = wb*(wh/wb)^((k+N+(1/2)*(1-r))/(2N+1))
%                 wk  = wb*(wh/wb)^((k+N+(1/2)*(1+r))/(2N+1))
%                 K   = wh^r.
    
    % Check number of input parameters
    if nargin < 4
        error('OUSTAFOD:NotEnoughInputArguments', ...
              'Not enough input arguments.');
    end
    
	mu=wh/wb;
	k=-N:N;
	w_kp=(mu).^((k+N+0.5-0.5*r)/(2*N+1))*wb;
	w_k=(mu).^((k+N+0.5+0.5*r)/(2*N+1))*wb;
	K=wh^r;
	G=tf(zpk(-w_kp',-w_k',K));
    
end


