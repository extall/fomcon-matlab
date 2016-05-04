function d2sos2(D)
%D2SOS Convert discrete form IIR filter into SOS form and echo the
%resulting coefficients in C style arrays to screen

[z,p,k] = zpkdata(D,'v');
[sos] = zp2sos(z,p,k);
 sosStr = sos2iir(sos);
 
 disp(char(13));
 
 disp(' ----- Coefficient arrays -----');
 disp(sosStr);
 
 disp(char(13));

end

