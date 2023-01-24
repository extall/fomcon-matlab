function display(p)
%DISPLAY Display the UFPOLY object
[str, uivl] = ufpoly2str(p.a, p.na, p.symb);
disp(' ');
disp(str);
disp(' ');
if uivl > 0
    disp('Fractional polynomial with uncertainty intervals.');
else
    disp('Fractional polynomial.');
end
disp(' ');
end

