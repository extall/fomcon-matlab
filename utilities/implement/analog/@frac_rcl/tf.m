function Z = tf(fract)
%TF Return Z(s) of the fractance circuit in form of a transfer function

    % Get function to handle the generation of the impedance
    str_fun = str2func(fract.structure);
    Z       = str_fun(fract.R, fract.C, fract.L, fract.K, fract.params);
    Z       = tf(Z);

end

