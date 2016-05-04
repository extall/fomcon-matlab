function fidata_plot_trim()
%FIDATA_PLOT_REMSOFF This is a helper function which facilitates removing
%                    the static offset in FIDATA structures output from the plot

    h      = gcf;
    useint = get(h, 'user');
    idd    = useint.idset;
    
    % Removing offset
    name='Trim range';
    
    prompt={'New workspace variable name:'};
    numlines=1;
    defaultAnswer = {'_rs'};
    options.Resize      = 'on';
    options.WindowStyle = 'normal';
    trimData=inputdlg(prompt,name,numlines,defaultAnswer,options);
    
    if ~isempty(trimData)
        
        id1 = idd;
		id1.y = id1.y-id1.y(1);
    
        % And store new system
        assignin('base', trimData{1}, id1);
    
        % Open new plot
        plot(id1);
        
    end

end

