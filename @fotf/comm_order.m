function ord = comm_order(G, type)
%COMM_ORDER Obtain commensurate order of fractional system.
%
% This function will compute the commensurate order of a fractional
% system in FOTF form with minimal resolution q=0.01 (default).

% Usage: ORD = COMM_ORDER(G, TYPE)
% TYPE can either be 'num' or 'den' to obtain commmensurate order of either
% numerator or denominator, if omitted, the function will return the system
% commensurate order.

    % Load configuration parameters
    config = fomcon('config');
    
    % Get the factor
    min_comm_order = config.Core.Commensurate_order.Min_comm_order;
    comm_factor = 1/min_comm_order;

    if nargin < 2

        a = G.na;
        b = G.nb;

        a1 = fix_s(comm_factor*a);
        b1 = fix_s(comm_factor*b);

        ab = [a1 b1];

        n = gcd(ab(1), ab(2));

        for i=3:length(ab)
            n = gcd(n, ab(i));
        end

        ord = n/comm_factor;
        
    elseif nargin == 2
        
        switch type
            
            case 'num'
                a=G.nb;
            case 'den'
                a=G.na;
            otherwise
                error('TYPE must be either ''num'' or ''den''.');
                
        end
        
        a1 = fix_s(comm_factor*a);
        
        % Check if array has a single element
        if numel(a1) == 1
            
            if a1(1) < 1, a1(1) = 1; end
            n = a1(1);
        
        else
            
            % Number of elements greater than 1
            n=gcd(a1(1),a1(2));

            for i=3:length(a1)
                n=gcd(n,a1(i));
            end
        end
        
        ord = n/comm_factor;
        
    end

end

