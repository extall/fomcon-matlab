function s1 = foss(G)
%FOSS Convert a commensurate-order fractional-order transfer function to state-space form.
%
% Usage: S1 = foss(TF1)
% where TF1 is a FOTF object.

    [num,den,q] = tfdata(G);
    [A,B,C,D] = tf2ss(num,den);
    s1 = foss(A,B,C,D,q);

end

