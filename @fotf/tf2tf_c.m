function f1 = tf2tf_c(G)
%TF2TF_C Convert a FOTF object into a CRONE FRAC_TF object
% 
% *** REQUIRES Object-Oriented CRONE toolbox ***

    % Convert directly
    f1 = frac_tf(frac_poly_exp(G.b,G.nb), frac_poly_exp(G.a,G.na));

end

