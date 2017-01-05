function c = cminus(a, b)
%CMINUS Minus operation for column vectors
%   C = CMINUS(A,B) returns the element-by-element difference vector
%                   of vectors A and B where both will be initially
%                   transposed if necessary

c = colv(a) - colv(b);

end

