function self = zscale(self, gains)
%SCALE Frequency scaling

    % Process all gains
    if iscell(gains)
        newgains = self.K;
        for n=1:length(gains)
            if ~isempty(gains{n})
                if ~isempty(self.R) && ~isempty(self.R{n})
                    self.R{n} = self.R{n}./gains{n};
                end
                if ~isempty(self.C) && ~isempty(self.C{n})
                    self.C{n} = self.C{n}.*gains{n};
                end
                if ~isempty(self.L) && ~isempty(self.L{n})
                    self.L{n} = self.L{n}./gains{n};
                end  
                newgains{n} = gains{n};
            end
        end
        self.K = newgains;
    else
        self.R = self.R./gains;
        self.C = self.C.*gains;
        self.L = self.L./gains;
        self.K = gains;
    end

end

