function fidata_plot_trim()
%FIDATA_PLOT_TRIM This is a helper function which facilitates trimming
%                 FIDATA structures from the plot

    h      = gcf;
    useint = get(h, 'user');
    idd    = useint.idset;
    
    % Trimming range
    name='Trim range';
    
    prompt={'New workspace variable name:', ...
            'Keep data from T1 [s]: ', ...
            'to T2 [s]:'};
    numlines=1;
    defaultAnswer = {'_t', '0', num2str(idd.t(end))};
    options.Resize      = 'on';
    options.WindowStyle = 'normal';
    trimData=inputdlg(prompt,name,numlines,defaultAnswer,options);
    
    if ~isempty(trimData)
        
        % Trim, if possible
        t1 = str2num(trimData{2});
        t2 = str2num(trimData{3});
        id1 = trim(idd, t1, t2);
    
        % And store new system
        assignin('base', trimData{1}, id1);
    
        % Open new plot
        plot(id1);
        
    end

end

