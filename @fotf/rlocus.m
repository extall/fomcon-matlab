function varargout = rlocus(G)
%RLOCUS Root locus of a fractional-order system
%   Usage: RLOCUS(G)    plots the root locus of a fractional-order system G

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Note: This function is adopted from J.T.Machado's numerical algorithm; %
%       Original copyright notice is presented below.                    %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% This program plots the root locus of a fractional order system
%
% The the width of the plot is defined by the variables: xmin, xmax, ymax
% The number of large cells for the numerical evaluation is variable: n
% The number of small cells within each large cell: ns
%
% The system is defined in syst.m
%
% 2011 J. Tenreiro Machado

% Input argument check
if nargin < 1
    error('RLOCUS:NotEnoughInputArguments', 'Not enough input arguments.');
end

% Load configuration parameters
config = fomcon('config');

xmin       = config.Core.Root_locus.X_min;               % Left limit
xmax       = config.Core.Root_locus.X_max;               % Right limit
ymax       = config.Core.Root_locus.Y_max;               % Upper limit
n          = config.Core.Root_locus.Coarse_grid_points;  % Large grid
ns         = config.Core.Root_locus.Fine_grid_points;    % Small grid

% Computation
ns1        = ns + 1;
ns2        = ns + ns1;
dx         =(xmax-xmin)/n;
dy         = 0.5*ymax/n;
dx2        = dx/2;
dxs        = dx/ns;
dy2        = dy/2;
dys        = dy/ns;
asarg_lim1 = 0.01;
asarg_lim2 = asarg_lim1*0.1;
ip         = 1;
x          = xmin;

while x<=xmax
    y=dy;
    while y<=ymax
        [mod_k,arg_k] = syst(G,x,y);
        asarg_k=abs(arg_k);
        if asarg_k<asarg_lim1
            asarg_min=asarg_k; arg_k_min=arg_k; asarg_k_min=asarg_k;
            
            % Small grid
            xsmin=x; ysmin=y;
            for i1=1:ns2
                for j1=1:ns2
                    k1(i1,j1)=0;
                end;
            end;
            i1min=ns1; j1min=ns1;
            k1(i1min,j1min)=1;
            flagk=0;
            while (xsmin<=x-dx2)||(xsmin>=x+dx2)|| ...
                    (ysmin<=y-dy2)||(ysmin>=y+dy2)|| ...
                    (flagk==0)
                for i1=1:3
                    xs(i1)=xsmin+(i1-2)*dxs;
                end;
                for j1=1:3
                    ys(j1)=ysmin+(j1-2)*dys;
                end;
                for i1=1:3
                    for j1=1:3
                        flagk=0;
                        i10=i1min+i1-2; j10=j1min+j1-2;
                        if (i10<=ns2) && (j10<=ns2) && (i10>0) && (j10>0)
                            if k1(i10,j10)==0
                                [mod_k,arg_k] = syst(G,xs(i1),ys(j1));
                                i1min=i10; j1min=j10;
                                k1(i10,j10)=1;
                                flagk=1;
                                asarg_k=abs (arg_k);
                                if asarg_k<asarg_min
                                    asarg_min=asarg_k; arg_k_min=arg_k;
                                    k_min=mod_k;
                                    xsmin=xs(i1); ysmin=ys(j1);
                                end;
                            end;
                        end;
                    end;
                end;
            end;
            
            % Small grid
            if abs(arg_k_min)<asarg_lim2
                ip     = ip+1;
                px(ip) = xsmin; py(ip) = ysmin;
                ip     = ip+1;
                px(ip) = xsmin; py(ip) = -ysmin;
            end;
        end;
        y=y+dy;
    end;
    x=x+dx;
end;

% Real axis
y=0.0;
x=xmin;
while x<=xmax
    [mod_k,arg_k] = syst(G,x,y);
    asarg_k=abs(arg_k);
    if asarg_k < asarg_lim2
        asarg_min = asarg_k; arg_k_min = arg_k; asarg_k_min = asarg_k;
        xsmin=x; ysmin=0;
        ip=ip+1;
        px(ip) = xsmin; py(ip) = ysmin;
    end;
    x=x+dx;
end;
if ip>1

    if nargout > 0
        Rts = [];
        Ks  = [];
    else
        h=gcf;
    end
    
    i=1;
    while i<=ip
        [qa,qna,qb,qnb] = fotfparam(G);
        root = px(i) + sqrt(-1)*py(i);
        k = 1/abs((qb*((root).^qnb.')/(qa*((root).^qna.'))));
        if nargout > 0
            Rts(end+1)=root; Ks(end+1)=k;
        else
            thisRoot = scatter(px(i), py(i), 3, [0 0 1], 'filled');
            params=struct; params.root = root; params.k=k;
            set(thisRoot, 'UserData', params);
            hold on;
        end
        i=i+1;
    end;
    
    if nargout > 0
        varargout{1} = Rts;
        varargout{2} = Ks;
    else
        % Set axes limits
        x_max = max(px); x_min = min(px);
        y_max = max(py); y_min = min(py);
        x_max = x_max + abs(x_max)/2; x_min = x_min - abs(x_min)/2;
        y_max = y_max + abs(y_max)/2; y_min = y_min - abs(y_min)/2;
        xlim([x_min x_max]); ylim([y_min y_max]);
        
        % Set axes labels
        xlabel ('Re(s)');
        ylabel ('Im(s)');
        grid on;
        hold off;
        datacursormode on;
        hh = datacursormode(h);
        set(hh,'UpdateFcn',@tooltip_text);
    end
end
end

% syst function prepares the data for calculation from a fotf object G
function [mod_k, arg_k] = syst(G,x,y)

% Numerator and denominator definitions
[b,bf,a,af] = fotfparam(G);
b=fliplr(b);bf=fliplr(bf);a=fliplr(a);af=fliplr(af);
imax=length(a); jmax=length(b);

% Total
logz=0.5*log(x*x+y*y);
argz=atan2(y,x);

% Numerator calculation
re (1)=0; im (1)=0;
for i=1:imax
    if af(i)==0, r=1; t=0; end
    if af(i)==1, r=exp(logz); t=argz; end
    if (af(i)~=0)&&(af(i)~=1), r=exp(af(i)*logz); t=af(i)*argz; end
    ar=a(i)*r;
    re(1)=re(1)+ar*cos(t); im(1)=im(1)+ar*sin(t);
end

m1=re(1)*re(1)+im(1)*im(1);

% Denominator calculation
re (2)=0; im (2)=0;
for j=1:jmax
    if bf(j)==0, r=1; t=0; end
    if bf(j)==1, r=exp (logz); t=argz; end
    if (bf(j)~=0)&&(bf(j)~=1), r=exp(bf(j)*logz); t=bf(j)*argz; end
    br=b(j)*r;
    re(2)=re(2)+br*cos(t); im(2)=im(2)+br*sin(t);
end

m2=re(2)*re(2)+im(2)*im(2);

% Output
mod_k=sqrt(m2/m1);
arg_k=atan2(re(1)*im(2)-im(1)*re(2),re(1)*re(2)+im(1)*im(2))+pi;
na=round(arg_k/6.283185307);
arg_k=arg_k-6.283185307*na;
end

% Add tooltips with relevant data to the plot
function txt = tooltip_text(empt, event_obj)
    % Load configuration parameters
    config = fomcon('config');
    numSigDig = config.Core.General.Model_significant_digits;
    
    pt   = get(event_obj, 'Target');
    params    = get(pt, 'UserData');
    txt = {['Root: ', num2str(params.root,numSigDig)], ...
           ['Gain: ', num2str(params.k,numSigDig)]};
end

