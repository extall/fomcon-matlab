function [G, J] = frid_prm(id1, Gini, ops)
%FRID_PRM Least squares freq. response fitting for fotf process models
% Use as follows: [G, J] = frid_prm(id1, Gini, ops)
% where G is the found tf, J is the final perf index (E of resids.^2)
% and id1 is the ffidata, Gini is a fotf in either 'fo' or 'so' form
% (see below), and ops are identification options.
%
%                         K
% 'fo' means G(s) = ------------ exp(-Ls)
%                      alpha 
%                   a s      + 1
%
%                                    K
% while 'so' means G(s) = ----------------------- exp(-Ls)
%                            beta       alpha 
%                         b s      + a s      + 1
%

% Check input system
[a,na,b,nb,L] = fotfparam(Gini);

if numel(b) == 1 && numel(a) == 2 && ~isempty(L)
    get_model = @(x) fotf([x(2) 1],[x(3) 0],x(1),0, x(4));
    x0 = [b(1) a(1) na(1) L];
elseif numel(b) == 1 && numel(a) == 3 && ~isempty(L)
    get_model = @(x) fotf([x(2) x(4) 1],[x(3) x(5) 0], x(1),0, x(6));
    x0 = [b(1) a(1) na(1) a(2) na(2) L];
else
    error('Unsupported model structure.');
end

% Get the magnitude
magn = id1.mag;
phas = id1.phase;

% Complex response and frequencies
r = 10.^(magn./20).*exp(deg2rad(phas).*1i);
w = id1.w;

% Define optimization problem
[x, J] = fminsearch(@(x) sumsq(fr_resid(x, get_model, r, w)), x0, ops);

% Return model
G = get_model(x);

end

function res = fr_resid(x, fun, r1, w)
    
    % Residuals are abs values + phase in deg
    r = squeeze(freqresp(fun(x), w)).';
    
    % Difference
    res = r - r1;
    
    % Distance between points in complex plane comprise the residuals
    for n=1:numel(res)
        res(n) = norm(res(n));
    end

end

function y = sumsq(x)
    y = sum(x.^2);
end