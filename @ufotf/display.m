function display(g)
%DISPLAY Display the UFOTF object

bp = g.b;
ap = g.a;
ioDelay = g.ioDelay;

[strb, uivlb] = ufpoly2str(bp.a, bp.na, 's');
[stra, uivla] = ufpoly2str(ap.a, ap.na, 's');

lensb = length(strb);
lensa = length(stra);

if lensb > lensa
    stra = make_samelen(stra, lensb);
elseif lensb < lensa
    strb = make_samelen(strb, lensa);
else
    % Do nothing, strings of equal length
end

sep = make_sep(max(lensa, lensb));

% Is there a delay term? If so, add it to sep
if ~isempty(ioDelay)
    delterm = num2str(ioDelay);
    if max(size(ioDelay)>1)
        delterm = ['-[' regexprep(delterm, '\s+', ', ') ']s'];
    else
        delterm = ['-' delterm 's'];
    end
    sep = [sep ' exp(' delterm ')'];
end

disp(' ');
disp(strb);
disp(sep);
disp(stra);
disp(' ');
if max(uivla, uivlb) > 0
    disp('Fractional-order transfer function with uncertainty intervals.');
else
    % TODO: Is this the optimal way?
    disp('Fractional-order transfer function member (convert to fotf)');
end
disp(' ');
end

function newstr = make_samelen(str, longest)
    n_spaces = floor((longest - length(str)) / 2);
    add_spaces = repmat(char(' '), 1, n_spaces);
    newstr = [add_spaces str add_spaces];
end

function sep = make_sep(longest)
    sep = repmat(char('-'), 1, longest);
end