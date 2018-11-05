function plot(idd)
%PLOT Plot the observed experimental results in the time domain.

    h = figure;
    set(h,'toolbar','figure');
    
    % Plot response
    ax1 = subplot(2,1,1);
    plot(idd.t, idd.y, 'Color', 'b', 'Linewidth', 2);
    ylabel('System output y(t)');
    grid;
    
    % Plot input signal
    ax2 = subplot(2,1,2);
    plot(idd.t, idd.u, 'Color', [0.25 0.25 0.25], 'Linewidth', 2);
    ylabel('System input u(t)');
    xlabel('Time [s]');
    grid;
    
    linkaxes([ax1 ax2], 'x');
    
    set(h, 'NumberTitle', 'off');
	set(h, 'Name', 'Identification data');
	
    % Create user interaction structure
    useint.idset = idd;
    set(h, 'user', useint);
    
    % Add menu
    f = uimenu('Label', 'Data', 'Position', 2);
        uimenu(f, 'Label', 'Trim', 'MenuSelectedFcn', 'fidata_plot_trim');
        uimenu(f, 'Label', 'Remove static offset', 'MenuSelectedFcn', ...
               'fidata_plot_remsoff');

end

