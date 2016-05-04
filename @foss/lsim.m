function varargout = lsim(S, u, t)
%LSIM Linear simulation of FOSS systems.
%
% Usage: Y = LSIM(S, U, T)

% Check input arguments
if nargin < 3
    error('LSIM:NotEnoughInputArguments', ...
        'Not enough input arguments.');
end

% Numbers of outputs and inputs
noOuts = size(S.fotf_array,1);
noIns  = size(S.fotf_array,2);

% Time vector length
tLen = numel(t);

if size(u,1) ~= tLen || size(u,2) ~= noIns
    error(['When simulating the response to a specific input ' ...
        'signal, the input data U must be a matrix with ' ...
        'as many rows as samples in the time vector T, and ' ...
        'as many columns as input channels.']);
end

y_cell = cell(noOuts, noIns);
y_array = zeros(tLen, noOuts);

% Go through all of the FO TFs in the FOSS array
for k=1:noOuts
    for l=1:noIns
        y_cell{k,l} = lsim(getfotf(S,k,l),u(:,l),t);
        y_array(:,k) = y_array(:,k) + y_cell{k,l};
    end
end

% Return the outputs array
varargout{1} = y_array;

% Also return individual FOTF outputs
if nargin > 1
    varargout{2} = y_cell;
end

end

