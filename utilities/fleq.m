function result = fleq(a,b)
%FLEQ Compare two floating-point numbers based on EPS
    result = abs(a-b) <= eps;
end

