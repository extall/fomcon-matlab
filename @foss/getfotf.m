function G = getfotf(S, k, l)
%GETFOTF Return FO TF to Kth output from Lth intput

if nargin < 1
    error('GETFOTF:NotEnoughInputArguments', ...
        'Not enough input arguments.');
end

if nargin < 2
    G = S.fotf_array{:};
elseif nargin < 3
    G = S.fotf_array{k, :};
elseif nargin >= 3
    G = S.fotf_array{k,l};
end

