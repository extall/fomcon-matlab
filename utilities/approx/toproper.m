function G = toproper(H, tau)
%TOPROPER Create proper LTI models by applying low-pass filters
%
%  Usage: G = TOPROPER(H, tau)
%  where  H   - a non-proper LTI system (with more zeros than poles),
%         tau - low-pass filter time constant.
%
%  Returns a ZPK object in form
%
%                                      n
%                    (      1        )
%         G = H(s) * ( ------------- ),
%                    ( 1 + (1/tau)*s )
%
%  where n is the number of added poles until the system is proper.


    % Check input
    if nargin < 2
        error('TOPROPER:NotEnoughInputArguments', 'Not enough input arguments');
    end
    
    % Convert system to ZPK, if needed
    if strcmpi(class(H),'tf'), H = zpk(H); end
    
    % Low-pass filter definition
    Hf = zpk(1/(1+tf('s')*(1/tau)));
    
    % Convert to state space?
    if strcmpi(class(H), 'ss'), Hf = ss(Hf); end
    
    % Make sure that the system is proper
    % (number of poles is greater than number of zeros)
    G = H;
    while ~isproper(G)
        G = G * Hf;
    end

end

