function varargout=lsim(G, u, t)
% LSIM Linear simulation of a fractional-order dynamic system.
%
% Usage: Y = LSIM(G, U, T) where G is the FOTF object to simulate,
%                          U is the input signal vector and T is the
%                          time sample vector

    varargout = {};

    % Check input arguments
    if nargin < 3
        error('LSIM:NotEnoughInputArguments', ...
		'Not enough input arguments.');
    end
    
    if min(size(u)) ~= 1 
        error('LSIM:UNotVector', 'u(t) is not a vector');
    end
    
    if min(size(t)) ~= 1
        error('LSIM:TNotVector', 't is not a vector');
    end

	% Determine progress bar requirement
	if varexists('pbar_show_') && evalin('base', 'pbar_show_')
		showPBar = true;
	else
		showPBar = false;
	end

    % Coefficient and exponent matrices
    a=G.a;
    eta=G.na;
    b=G.b;
    gamma=G.nb;

    nA=length(a);
    h=t(2)-t(1);
    D=sum(a./[h.^eta]);
    W=[];
    nT=length(t);
    vec=[eta gamma];
    D1=b(:)./h.^gamma(:);
    y1=zeros(nT,1);
    W=ones(nT,length(vec));
    
    if showPBar, wBar = pbar('Computation step 1 of 3...'); end
    
    for j=2:nT
        W(j,:)=W(j-1,:).*(1-(vec+1)/(j-1));
        if showPBar, update(wBar, j, nT); end
    end

    if showPBar, wBar = pbar('Computation step 2 of 3...'); end
    
    for i=2:nT
        A=(y1(i-1:-1:1))'*W(2:i,1:nA);
        y1(i)=(u(i)-sum(A.*a./(h.^eta)))/D;
        if showPBar, update(wBar, i, nT); end
    end
    
    if showPBar, wBar = pbar('Computation step 3 of 3...'); end
    
    y = zeros(nT,1);    % Initialize y
    for i=2:nT
        y(i)=(W(1:i,nA+1:end)*D1)'*(y1(i:-1:1));
        if showPBar, update(wBar, i, nT); end
    end
	
    % Account for I/O delay
    if (~isempty(G.ioDelay) && G.ioDelay > 0)
        ii = find(t>G.ioDelay);
        
        % There is a possibility that the value of the sampling interval
        % is greater or equal to the delay. In this case we disregard it.
        if ~isempty(ii)
            lz = zeros(ii(1)-1,1);
            y  = [lz; y(1:end-length(lz))];
        end
    end
    
	if showPBar, delete(wBar); end
    
    % Make sure y is a column vector
    if size(y, 1) < size(y, 2), y=y'; end
    
    switch (nargout)
    case 0
        plot(t,y);
        title('Linear system simulation');
        xlabel('Time [s]');
        ylabel('Amplitude');
    case 1
        varargout{1} = y;
    end
    
end