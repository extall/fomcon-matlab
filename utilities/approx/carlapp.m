function H = carlapp(M, N, G, G0, retol)
%CARLAPP Carlson's approximation of the fractional-order transfer function G(s)^(1/M)
%
% Usage: H = CARLAPP(M, N, G, G0, RETOL),
%
%      where M is the denominator of ALPHA = 1/M, where ALPHA is
%            approximated by G(s)^ALPHA and N is the order of
%            approximation,
%			 G is the optional function G(s) taken as G(s) = s by default,
%            G0 is the initial solution, if omitted or empty considered to
%            be zpk('1'),
%            RETOL is the tolerance for model reduction using a matching
%            zero-pole reduction method realized by applying the minreal()
%            function.

	% Check input
    if nargin < 1
		error('CARLAPP:NotEnoguhInputArguments', ...
              'Not enough input arguments.');
    end
    
    % Check approximation order
    if nargin < 2
        N = 2;
    end
	
    % Check transfer function
	if nargin < 3 || isempty(G), G = zpk('s'); end
    
    % Check initial solution
    if nargin < 4 || isempty(G0), G0 = zpk('1'); end
    
    % Check model reduction tolerance
    if nargin < 5
        retol = [];
    end

    % Ensure system is ZPK
	G = zpk(G);
	
	% Initial transfer function
	H = zpk(G0);
    
    % Begin iterations
    for n=1:N
        
        % Carlson formula
		H = H * ((M+1)*G + (M-1)*H^M) / ...
                ((M-1)*G + (M+1)*H^M);
            
        % Model reduction
        if ~isempty(retol)
            H = minreal(H, retol);
        end
        
    end

end



