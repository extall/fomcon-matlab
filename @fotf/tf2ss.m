function [A,B,C,D,q] = tf2ss(G)
%TF2SS Convert a commensurate-order fractional-order transfer function to state-space form.
%
% Usage: [A,B,C,D,q] = tf2ss(TF1)
% where TF1 is a FOTF object.

    [num,den,q] = tfdata(G);
    [A,B,C,D] = tf2ss(num,den);
	
end

