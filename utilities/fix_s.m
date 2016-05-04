function a = fix_s(b)
%FIX_S Extract integer part from a given number by using sprintf
    for n=1:length(b)
        a(n) = fix(str2num(sprintf('%d',b(n))));    
    end
end