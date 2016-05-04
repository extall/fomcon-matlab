function [a, na, b, nb, G, J] = ffid(idata, q, ord, method)
%FFID Frequency-domain identification of fractional-order models
%
% This function uses methods from NINTEGER toolbox by Duarte Valerio.
%
% Usage: [A, NA, B, NB, G, J] = FFID(IDATA, Q, ORD, METHOD)
%
% where   A, NA, B, NB, G, J - identified model coefficients/exponents,
%                               corresponding FOTF object and J error index
%
%         IDATA - FFIDATA structure with observed magnitude, phase and 
%                 frequencies at which the observations were made
%         Q     - commensurate order of model
%         ORD   - orders of pole and zero polynomials: [N; M]; for
%                 'hartley' method only the pole polynomial N is required.
%         method- identification method: 'h' - Hartley, 'l' - Levy or 'v' -
%                 Vinagre (default). Please consult the NINTEGER toolbox
%                 documentation for details about each method.
    
    % Check input arguments
    if nargin < 4
        method = 'v';
    end
    
    if nargin < 3
        error('FFID:NotEnoughInputArguments', 'Not enough input arguments.');
    end
    
    % Get orders
    if strcmpi(method, 'h')
        n = ord(1);
    else
        n = ord(1);
        m = ord(2);
    end
    
    % Check commensurate order
    if q < 0 || q > 2
        error('FFID:InvalidCommensurateOrder', 'Commensurate order q must be in range 0<q<2');
    end
    
    % Identify model
    switch(lower(method))
        case 'h'
            g = hartley(idata.w, idata.mag, idata.phase, q, n, 'den');
            
            % Zero polynomial
            b = 1;
            nb = 0;
            
            % Pole polynomial
            a = real(g)';
            na = q*(n:-1:0);
            
        case 'l'
            g = levy(idata.w, idata.mag, idata.phase, q, n, m);
            
            b = g.num';
            nb = q*(m:-1:0);
            a = g.den';
            na = q*(n:-1:0);
           
        case 'v'
            g = vinagre(idata.w, idata.mag, idata.phase, q, n, m);
            
            b = g.num';
            nb = q*(m:-1:0);
            a = g.den';
            na = q*(n:-1:0);
    end
    
    % Build model
    G = fotf(a, na, b, nb);
    
    % Calculate error index
    if nargout > 5
        
        if strcmpi(method, 'h')
            tnum = 1;
            g1.den = g;
        else
            tnum = 0;
            for c = 1:length(g.num)
                tnum = tnum + g.num(c) * (1i*idata.w).^((m-c+1)*q);
            end
            g1.den = g.den;
        end
        
        tden = 0;
        for c = 1:length(g1.den)
            tden = tden + g1.den(c) * (1i*idata.w).^((n-c+1)*q);
        end
        
        J = sum((abs(10.^(idata.mag/20).*exp(1i*deg2rad(idata.phase)) - tnum./tden)).^2) / length(idata.w);
        
    end

end

