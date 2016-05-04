function s1 = tf2ss_c(G)
%TF2SS_C Convert a commensurate-order fractional-order transfer function to CRONE state space.
%
% Usage: S1 = tf2ss_c(TF1)
% where TF1 is a FOTF object
%
% *** REQUIRES Object-Oriented CRONE toolbox ***

    tf1 = frac_tf(frac_poly_exp(G.b,G.nb),frac_poly_exp(G.a,G.na));
    s1 = tf2ss(tf1);

end

