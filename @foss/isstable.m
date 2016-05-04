function K = isstable(S, doPlot)
%ISSTABLE Check stability of FOSS system.

% Check input
if nargin < 2
    doPlot = false;
end

q = S.q;
p = eig(S.A);
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

