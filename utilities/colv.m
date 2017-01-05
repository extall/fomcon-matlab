function a = colv(b)
%COLV Convert to column vector
%   A = COLV(B) converts B to a column vector assuming B is a vector

% Check whether input is a vector
[yb,xb] = size(b);
if min([yb,xb]) > 1
    error('COLV:OnlyVectorsSupported', ...
        'Only vectors are supported by this function.');
end

% Initial assignment
a = b;

% Transpose if row vector
if xb > yb
    a = a.';
end

end

