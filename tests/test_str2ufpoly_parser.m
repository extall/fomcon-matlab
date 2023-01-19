%% Testcase
mystr = '-3s^[1.3,1.5] + [12; 15] s^1.3e-9 - [9, 11]';
[a, na, symb] = str2ufpoly(mystr);

assert(all(a == [-3, -3; 12 15; -11, -9], 'all'), ...
    'Incorrect values of extracted coefficients');

assert(all(na == [1.3, 1.5; 1.3e-9 1.3e-9; 0, 0], 'all'), ...
    'Incorrect values of extracted exponents');

assert(symb == 's', ...
    'Incorrect extracted symbol');

%% Testcase
mystr = 'x+[2;3]';
[a, na, symb] = str2ufpoly(mystr);

assert(all(a == [1, 1; 2 3], 'all'), ...
    'Incorrect values of extracted coefficients');

assert(all(na == [1, 1; 0, 0], 'all'), ...
    'Incorrect values of extracted exponents');

assert(symb == 'x', ...
    'Incorrect extracted symbol');

%% Testcase
mystr = '[-1.1, +1.15]x+[2;3]';
[a, na, symb] = str2ufpoly(mystr);

assert(all(a == [-1.1, 1.15; 2 3], 'all'), ...
    'Incorrect values of extracted coefficients');

assert(all(na == [1, 1; 0, 0], 'all'), ...
    'Incorrect values of extracted exponents');

assert(symb == 'x', ...
    'Incorrect extracted symbol');

%% Testcase
mystr = '- [9, 11]';
[a, na, symb] = str2ufpoly(mystr);

assert(all(a == [-11, -9], 'all'), ...
    'Incorrect values of extracted coefficients');

assert(all(na == [0, 0], 'all'), ...
    'Incorrect values of extracted exponents');

assert(symb == 's', ...
    'Incorrect extracted symbol');