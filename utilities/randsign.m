function r = randsign(varargin)
%RANDSIGN Random sign generator
% The syntax of this function is the same as for the rand() function,
% however it returns vectors/matrices/arrays of random signs (+/-1).
%
% See also: rand

r = rand(varargin{:});
r(r<0.5)=-1; r(r>0.5)=1;

end

