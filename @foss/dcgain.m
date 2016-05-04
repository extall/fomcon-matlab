function K = dcgain(S)
%DCGAIN Static gain of the FO SS
%   K = dcgain(S) computes the steady-state gains for FOSS.

% Get necessary array size
noOuts = size(S.fotf_array,1);
noIns  = size(S.fotf_array,2);

% Create the response array
K = zeros(noOuts,noIns);

for k=1:noOuts
    for l=1:noIns
        K(k,l) = dcgain(getfotf(S,k,l));
    end
end

end

