function [qv, expr, tex] = polycfe(varargin)
%POLYCFE Polynomial continued fraction expansion
%   Usage: R = POLYCFE(B, A)
%          R = POLYCFE(G)
    
	config = fomcon('config');
	epsi   = config.Core.General.Internal_computation_accuracy;
	
    % Check input arguments
    if nargin == 1
       [b, a] = tfdata(varargin{1}, 'v');
    elseif nargin >= 2
       b = varargin{1};
       a = varargin{2};
    end

    % Tolerance
    eps2 = epsi;
    % Remove leading zeros
    b=b(find(abs(b)>eps2,1,'first'):end);
    a=a(find(abs(a)>eps2,1,'first'):end);

    % Continued fraction coefficients in s^n, n=0, 1, ...
    qv = {}; qq = {};

    % This will only be used if polynomials need to be inverted initially
    c = [];
    
    % Check degrees
    if length(b)<length(a)
       c = b; b = a; a = c;
       warning('POLYCFE:SwappingNumDen', ...
           'Swapping polynomial numerator and denominator.');
    end

    while numel(a) > 0

        % Polynomial long division
        [q, r] = deconv(b, a);

        % Take highest degree term in q
        hi  = q(1);
        deg = length(q)-1;

        % If this is not the only term, perform additional
        % subtraction to arrive at the necessary remainder 
        if length(q) > 1
            r1 = conv([hi zeros(1,deg)], a);
            r = leftpadz(b,r1) - leftpadz(r1,b);
        end

        % Add transfer function to list(s) only if the vector is nonzero
        if max(abs(hi)>eps2) > 0
            elem = zeros(1, deg+1);
            elem(1) = hi;
            qv{end+1} = elem;
        end
        
        if deg>0, deg_str = ' s'; else deg_str = ''; end
        if deg>1, deg_str = [deg_str '^' num2str(deg)]; end
        qq{end+1} = [num2str(hi) deg_str];

        % Swap polynomials
        b = a;
        a = r;

        % Remove leading zeros
        b=b(find(abs(b)>eps2,1,'first'):end);
        a=a(find(abs(a)>eps2,1,'first'):end);

    end
    
    % Form expression
    if nargout > 1
       
        expr_string=['1/(' strrep(qq{end},' ','*') ')'];
        
        for n=(length(qq)-1):-1:1
            part_expr = [strrep(qq{n},' ','*') '+' expr_string];
            if n~= 1
               expr_string = ['1/(' part_expr ')'];
            else
               expr_string = part_expr; 
            end
        end
        
        if ~isempty(c)
            expr_string = ['1/(' expr_string ')'];
        end
        
        expr = expr_string;
        
    end
    
    % Form TeX string using \cfrac
    if nargout > 2
        
        % Form \cfrac string
        tex_string = ['\cfrac{1}{' strrep(qq{end},' ','*') '}'];
        
        for n=(length(qq)-1):-1:1
            part_string = [strrep(qq{n},' ','') '+' tex_string];
            if n~= 1
                tex_string = ['\cfrac{1}{' part_string '}'];
            else
                tex_string = part_string;
            end
        end
        
        % Need to start with fraction
        if ~isempty(c)
            tex_string = ['\cfrac{1}{' tex_string '}'];
        end
        
        tex = tex_string;
        
    end
    
    % First term is zero
    if ~isempty(c), qq = {'0' qq{:}}; end
        
end

% Helper function for polynomial addition
function v = leftpadz(p1, p2)
    v = [zeros(1,max(0,numel(p2) - numel(p1))),p1];
end


