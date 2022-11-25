function varargout = fotf2expr(G)
%FOTF2EXPR converts a FOTF object to a symbolic expression
%   If one output argument is provided, returns a single expression.
%   If two or three, then returns the numerator, denominator, and ioDelay
%   as separate terms. If NARGIN=2 and ioDelay!=0, ioDelay will be included
%   in the numerator expression.

syms('s');
a = G.a; na = G.na;
b = G.b; nb = G.nb;
ioDelay = G.ioDelay;

% Create the expressions in string form first
sa  = fpoly2str(a, na, 's', '*');
sb  = fpoly2str(b, nb, 's', '*');
sio = iod2str(ioDelay, 's', '*');

% Replace empty delay with a mul by 1
if isempty(sio)
    sio = '1';
end

% Create actual symbolic expressions
expra   = eval(['(' sa ')']);
exprb   = eval(['(' sb ')']);
expriod = eval(sio);

% Depending on how many output arguments were requested, return the items
if nargout == 1
    varargout{1} = expra/exprb * expriod;
elseif nargout == 2
    varargout{1} = exprb * expriod;
    varargout{2} = expra;
elseif nargout == 3
    varargout{1} = exprb;
    varargout{2} = expra;
    varargout{3} = expriod;
end

end

