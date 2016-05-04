function [fractf] = oustapp(G, wb, wh, N, method)
%OUSTAPP Obtain a integer-order approximation of a fractional-order system.
%
% Usage: TF = OUSTAPP(G, wb, wh, N, method) 
%        
%        where G - fotf object,
%              wb, wh - lower and higher frequency bounds
%              N - approximation order
%              method - 'oust' for Oustaloup's method (default) or 
%                       'ref' for the refined Oustaloup method
%
%              (defaults: wb = 0.001, wh = 1000, N=5, method='oust')
%
%        See also: fsparam

    % Load configuration parameters
    config = fomcon('config');
    ft_d   = config.Approximations.Oustaloup_filter.Default_filter_type;
    wb_d   = config.Approximations.Oustaloup_filter.Default_low_freq_bound;
    wh_d   = config.Approximations.Oustaloup_filter.Default_high_freq_bound;
    N_d    = config.Approximations.Oustaloup_filter.Default_approx_order;

    if nargin < 5
        method=ft_d;
    end
    
    % Order of approximation
    if nargin < 4
        N=N_d;
    end

    % Higher frequency bound
    if nargin < 3
        wh = wh_d;
    end

    % Lower frequency bound
    if nargin < 2
        wb = wb_d; 
    end
    
    if nargin < 1
        error('OUSTAPP:NotEnoughInputArguments', ...
              'Not enough input arguments.');
    end
    
    method = lower(method);

    switch method
        case 'oust'
            doApp = 0;
        case 'ref'
            doApp = 1;
        otherwise
            warning('OUSTAPP:BadMethod', ...
                    'Method must be ''oust'' or ''ref'', using ''oust'' as default');
            doApp = 0;
    end
    
    s = tf('s');
    
    % Get FOTF object parameters
    [a,na,b,nb,ioDel] = fotfparam(G);
    
    % Create zero and pole arrays
    zeroPoly = 0;
    polePoly = 0;
    
    % Go through zero array
    for n=1:length(b)
       
        thisExp = nb(n);     % Current exponent
        
        % DEBUG
        intPart = fix_s(thisExp);           % Integer part of the exponent
        nintPart = double(thisExp-intPart); % Fractional part of the exponent
        
        toAdd = (s^intPart)*b(n);
        
        if nintPart ~= 0
            
            if doApp == 1
                toAdd = toAdd * tf(new_fod(nintPart,N,wb,wh));
            else
                toAdd = toAdd * oustafod(nintPart,N,wb,wh);
            end
            
        end
        
        zeroPoly = zeroPoly + toAdd;
      
    end
    
    % Go through pole array
    for n=1:length(a)
       
        thisExp = na(n);     % Current exponent
        
        intPart = fix_s(thisExp);           % Integer part of the exponent
        nintPart = double(thisExp-intPart); % Fractional part of the exponent
        
        toAdd = (s^intPart)*a(n);
        
        if nintPart ~= 0
            if doApp == 1
                toAdd = toAdd * tf(new_fod(nintPart,N,wb,wh));
            else
                toAdd = toAdd * oustafod(nintPart,N,wb,wh);
            end
        end
        
        polePoly = polePoly + toAdd;
      
    end
    
    % Convert to ZPK model
    zeroPoly = zpk(zeroPoly);
    polePoly = zpk(polePoly);
    
    fractf = zeroPoly/polePoly;
    
    % Add ioDelay if present
    if ioDel > 0
       fractf.ioDelay = ioDel; 
    end
    
    if isempty(fractf.z{1}) && ...
            isempty(fractf.p{1}) && ...
            isnan(fractf.k)
        % Resulting model is wrong
        error('OUSTAPP:BadModel', ...
            ['Model is erroneous! Conversion aborted; ' ...
             'possibly due to too high order of approximation']);
    end
    
end

