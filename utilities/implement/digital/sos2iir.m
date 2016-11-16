function str = sos2iir(s1)
%SOS2IIR Creates a C style array for the IIR filter in SOS form
    
    str1 = '';
    str2 = '';
    
    for n=1:size(s1,1)
        str1 = [str1 '{' sprintf('%+.10f', s1(n,1)) ', ' ...
            sprintf('%+.10f', s1(n,2)) ', ' sprintf('%+.10f', s1(n,3)) '}'];
        if (n ~= size(s1,1)) str1 = [str1 ',' char(13)]; end
        % We do NOT show the +1 coefficient for discrete poles here
        str2 = [str2 '{' sprintf('%+.10f', s1(n,5)) ', ' ...
            sprintf('%+.10f', s1(n,6)) '}'];
        if (n ~= size(s1,1)) str2 = [str2 ',' char(13)]; end
    end
    
    str = [str1 char(13) char(13) str2];

end

