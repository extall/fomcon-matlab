function h = plot(st)
%PLOT Plots the generated signal sequence
h = figure; plot(st.t, st.u); xlabel('Time [s]'); ylabel('Amplitude');

end

