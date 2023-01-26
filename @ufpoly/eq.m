function is_equal = eq(p1, p2)
%EQ Check if the two objects are equal

% If sizes of matrices a1 and a2 are of different size,
% then clearly the objects are not equal.
if ~all(size(p1.a) == size(p2.a))
    is_equal = 0;
else
    is_equal = all(p1.a == p2.a, 'all') && ...
        all(p1.na == p2.na, 'all') && ...
        p1.symb == p2.symb;
end

