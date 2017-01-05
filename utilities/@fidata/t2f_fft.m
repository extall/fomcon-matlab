function id2 = t2f_fft(id1)
%T2F_FFT Convert time domain fidata to freq. domain ffidata with user input
%   The conversion is not yet fully automated. The user must locate the
%   meaningful frequency range manually

[yr, ur, w] = get_fft_fidata(id1);

r = yr./ur;
idt = ffidata(r,w);
h = figure;
bode(idt);

% The input
prompt = {'Start freq. range [rad/s]', 'End freq. range [rad/s]'};
dlg_title = 'Locate frequency range containing meaningful information.';
num_lines = 1;
defaultans = {'0.01', '100'};
ops = struct;
ops.WindowStyle = 'normal';
answer = inputdlg(prompt, dlg_title, num_lines, defaultans, ops);

if ~isempty(answer)
    w_lo = str2double(answer{1});
    w_hi = str2double(answer{2});
    r(w<w_lo) = []; w(w<w_lo) = [];
    r(w>w_hi) = []; w(w>w_hi) = []; 
end

% Return the corrected data set and close the window
if ishandle(h), delete(h); end
id2 = ffidata(r,w);

end

