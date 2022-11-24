function varargout=lsim(G, u, t, opts)
% LSIM Linear simulation of a fractional-order dynamic system.
%
% Usage: Y = LSIM(G, U, T, OPS) where
%    G is the FOTF object to simulate,
%    U is the input signal vector
%    T is the time sample vector
%    OPTS is a struct with configurable options:
%      OPTS.GL_Order = 1..3 is the GL simulation order (default: 1)
% 

    varargout = {};

    % Check input arguments
    if nargin < 3
        error('LSIM:NotEnoughInputArguments', ...
		'Not enough input arguments.');
    end
    
    options = []; % Empty matrix by default
    if nargin > 3
        if ~isstruct(opts)
            error('LSIM:WrongOptionsArgument', ...
                'The OPTS argument must be a struct containing config options.');
        else
            options = opts; % Now a structure
        end
    end
    
    if min(size(u)) ~= 1 
        error('LSIM:UNotVector', 'u(t) is not a vector');
    end
    
    if min(size(t)) ~= 1
        error('LSIM:TNotVector', 't is not a vector');
    end
    
    % For each configurable parameter, check the options
    conf = fomcon('config');
    
    % Simulation precision. Depending on this parameter, the specific
    % algorithm used for computations will also be different.
    GL_Order = 1; % Default solver
    if cfieldexists(options, 'GL_Order')
        GL_Order = options.GL_Order;
    else
        GL_Order = conf.Core.Simulations.GL_Order;
    end
    
    if ~isnumeric(GL_Order) || GL_Order < 1 || GL_Order > 3
        error('LSIM:SimOrderOutOfRange', ...
            'The GL_Order parameter must be in the range 1..3');
    end
    
	% Determine progress bar requirement
    % TODO: display of progress bar should be a global option instead
    % TODO: add this as a config parameter
    if varexists('pbar_show_') && evalin('base', 'pbar_show_')
		showPBar = true;
    else
		showPBar = false;
    end
    
    % Decide which algorithm to use
    if GL_Order == 1

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
    
    else
        
        % Extract the information from the fotf object
        a = G.a;
        na = G.na;
        b = G.b;
        nb = G.nb;
       
        % Use the more precise solver adopted from (Xue 2017)
        h = t(2)-t(1);
        n = length(t);
        vec = [na nb];
        u = colv(u); t = colv(t);
        g = gene_params(GL_Order);
        W = [];
        
        % Get the weighting coefficients
        if showPBar, wBar = pbar('Computation step 1 of 3...'); end
        lv = length(vec);
        for i=1:lv
            W = [W; get_vecw(vec(i), n, g)];
            if showPBar, update(wBar, i, lv); end
        end
        
        D1 = colv(b)./h.^colv(nb);
        nA = length(a);
        y1 = zeros(n,1);
        W = W.';
        D = sum((a.*W(1,1:nA))./(h.^na));
        
        % Reaction to input
        if showPBar, wBar = pbar('Computation step 2 of 3...'); end
        for i=2:n
            A=[y1(i-1:-1:1)].'*W(2:i,1:nA);
            y1(i)=(u(i)-sum(A.*a./(h.^na)))/D;
            if showPBar, update(wBar, i, n); end
        end
        
        % Final computations
        if showPBar, wBar = pbar('Computation step 3 of 3...'); end
        for i=2:n
            y(i) = (W(1:i, nA+1:end) * D1).'*(y1(i:-1:1));
            if showPBar, update(wBar, i, n); end
        end
        
    end

    % Make sure y is a column vector
    y = y(:);
	
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

% Generating function coefficients necessary to do numerical
% simulations with errors o(h^2), o(h^3) and so on.
function params = gene_params(order)
    G = { ...
      [3/2 -2 1/2], ...
      [11/6 -3 3/2 -1/3], ...
    }; % Can be continued, if needed; or use general result from Xue2017.
    params = G{order-1};
end

% Weighting coefficient computation for higher order simulations
function w=get_vecw(gam, n, g)
% Copyright (c) Dingyu Xue, Northeastern University, China
%
% For further clarifications, please see the FOTF toolbox where this
% function was taken from:
% https://se.mathworks.com/matlabcentral/fileexchange/60874-fotf-toolbox

p=length(g)-1; b=1+gam; g0=g(1);
w=zeros(1,n); w(1)=g(1)^gam;
for m=2:p, M=m-1; dA=b/M;
   w(m)=-[g(2:m).*((1-dA):-dA:(1-b))]*w(M:-1:1).'/g0;
end
for k=p+1:n, M=k-1; dA=b/M;
   w(k)=-[g(2:(p+1)).*((1-dA):-dA:(1-p*dA))]*w(M:-1:(k-p)).'/g0;
end
end
