function [K,q,err,apol]=isstable(G, doPlot)
% ISSTABLE Check fractional-order system stability

% Usage: [K,q,err,apol]=isstable(G, doPlot)
%
%        K - boolean, true if system is stable, false if not
%        q - calculated commensurate-order
%        err - stability assessment error norm
%        apol - closest pole absolute imaginary part value
%               to the unstable region
%
%        doPlot : if set to true, will create a plot with relevant system
%        poles

    % Load configuration parameters
    config = fomcon('config');
    epsi   = config.Core.General.Internal_computation_accuracy;
	
    % Get the factor
    min_comm_order = config.Core.Commensurate_order.Min_comm_order;
    comm_factor = 1/min_comm_order;

    % Check input
    if nargin < 2
        doPlot = false;
    end

    % No support for systems with delays :: TODO
    %if G.ioDelay ~= 0
    %    error('ISSTABLE:SystemWithDelaysNotSupported', ...
    %          'Cannot assess stability of a system with internal delays.');
    %end
    
    a=G.na;
    a1=fix_s(a*comm_factor);
    q = comm_order(G,'den');
    a=fix_s(a1/(comm_factor*q));
    b=G.a;
    
    c(a+1)=b;
    c=c(end:-1:1);
    p=roots(c);
    
    %Check if there is a single root
    if p == 0
        % Do nothing
    else
        p=p(abs(p)>epsi);
    end
    
    err=norm(polyval(c,p));
    apol=min(abs(angle(p)));
    K=apol>q*pi/2;
    
    % Check if drawing is requested
    if doPlot
        
        % Create new figure and plot poles
        h=gcf();
        plot(real(p),imag(p),'x',0,0,'o')
        
        % Get x axis limit
        xm=xlim;
        xm(1)=0;

        qpi = q*pi/2;

        x_fill = [xm(1) xm(2) xm(2) xm(1)];
        y_fill = [0 qpi -qpi 0];

        % Draw unstable region
        patch(x_fill, y_fill, 'r', 'FaceAlpha', 0.5);
        
    end
    
end