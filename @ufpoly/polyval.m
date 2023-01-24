function y = polyval(p, x)
%POLYVAL Compute the values of a member of a ufpoly

[str, unval] = ufpoly2str(p.a, p.na, p.symb, '*');

% Also check if a specific member is provided
if unval > 1
    error('Please provide a specific member of a UFPOLY using sample(p).');
end

% Evaluate to get the complete expression
syms(p.symb);

% Get a MATLAB function
mfunc = matlabFunction(eval(str));

% Finally, compute the values
y = mfunc(x);

end

