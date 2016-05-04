function Z = zpk(fract)
%ZPK Return Z(s) of the fractance circuit in form of a zero-pole-gain model

    % Get function to handle the generation of the impedance
    str_fun = str2func(fract.structure);
    Z       = str_fun(fract.R, fract.C, fract.L, fract.K, fract.params);
    Z       = zpk(Z);

end