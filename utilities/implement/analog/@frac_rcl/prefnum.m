function self = prefnum(self, Rp, Cp, Lp, Kp)
%PREFNUM Match components of the network to a set preferred numbers

    % Substitute preferred numbers depending on the number
    % of input arguments and the existence of these
    if nargin >= 2
        if ~isempty(Rp) && ~isempty(self.R)
            self.R = prefnumerize(self.R, Rp);
        end
    end
    
    if nargin >= 3
        if ~isempty(Cp) && ~isempty(self.C)
            self.C = prefnumerize(self.C, Cp);
        end
    end
    
    if nargin == 4
        if ~isempty(Lp) && ~isempty(self.L)
            self.L = prefnumerize(self.L, Lp);
        end
    end
    
    if nargin == 5
        if ~isempty(Kp) && ~isempty(self.K)
            self.K = 1/prefnumerize(1/self.K, Kp);
        end
    end
            
end

