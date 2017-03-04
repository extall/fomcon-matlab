function a = colv(varargin)
%COLV Convert to column vector
%   A = COLV(B,C,D,...) converts input vectors to a single column vector

a = [];

% Check whether input is a vector
for k=1:numel(varargin)
    
    b = varargin{k};
    [yb,xb] = size(b);
    if min([yb,xb]) > 1
        error('COLV:OnlyVectorsSupported', ...
            'Only vectors are supported by this function.');
    end
    
    % Transpose if row vector
    if xb > yb
        b = b.';
    end
    
    % Add to final vector
    a = [a; b];
end

end

