function Z = doustafod(r,N,wb,wh,Ts)
%DOUSTAFOD Discrete Oustaloup approximation of a s^r operator
%
% USAGE: Z = DOUSTAFOD(R,N,WB,WH,TS)
% where  R is the operator order,
%        N relates to the number of zeros and poles: N_ZP = 2N+1,
%        WB and WH are lower and higher frequency bounds of approximation,
%        TS is the sampling time.

    % Check number of input arguments
    if nargin < 5
        error('DOUSTAFOD:NotEnoughInputArguments', ...
              'Not enough input arguments.');
    end
    
    % Check upper frequency bound
    if (wh > 2/Ts)
        wh = 2/Ts;
        warning('DOUSTAFOD:UpperFrequencyBoundChanged', ...
                'The upper frequency bound has been changed to 2/Ts.');
    end
    
    % Compute continuous zeros/poles
    mu   =  wh/wb;
    k    = -N:N;
    w_kp = (mu).^((k+N+0.5-0.5*r)/(2*N+1))*wb;
    w_k  = (mu).^((k+N+0.5+0.5*r)/(2*N+1))*wb;

    % Discrete-time mapping
    zz = exp(-Ts*w_kp);
    zp = exp(-Ts*w_k);
    
    % Correction frequency and new system gain
    wu  = sqrt(wb*wh);
    Kdc = abs(prod(exp(sqrt(-1)*wu*Ts)-zz)/prod(exp(sqrt(-1)*wu*Ts)-zp));
    
    % Correct gain and correcting factor
    Ks = wu^r;
    Kc = Ks / Kdc;
    
    % Return ZPK object
    Z = zpk(zz,zp,Kc,Ts);