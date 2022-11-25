function varargout = fotf_gui(varargin)
% FOTF_GUI Fractional-order transfer function analysis tool
%
% The FOTF_GUI (alias FOMCON) tool allows to work with fractional-order
% transfer function models and serves as the starting point for associated
% workflows.

% Last Modified by GUIDE v2.5 03-Oct-2013 16:43:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fotf_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @fotf_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before fotf_gui is made visible.
function fotf_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fotf_gui (see VARARGIN)

% Choose default command line output for fotf_gui
handles.output = hObject;

% Initialize listbox
refreshFotfList(handles);

% Update handles structure
guidata(hObject, handles);

% Create figures storage
guiData.Fig_time_response = [];
guiData.Fig_freq_response = [];
set(hObject, 'UserData', guiData);

% UIWAIT makes fotf_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fotf_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnStable.
function btnStable_Callback(hObject, eventdata, handles)
% hObject    handle to btnStable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    selected = get(handles.fotfList, 'Value');
    allEntries = get(handles.fotfList, 'String');
    thisEntry = allEntries{selected};
    
    [varex, varcl] = varexists(thisEntry);
    
    if varex && varcl
        
        myFotf = evalin('base', thisEntry);
        
        % Create figure
        h  = figure;
        ax = axes;
        
        % Check stability
        disp(['Stability test for ''' thisEntry '''']);
        [K,q,err,apol] = isstable(myFotf, true)
        
		config = fomcon('config');
        numSigDig = config.Core.General.Model_significant_digits;
		
        % Display message box with resolution
        if (K == 1)
            msgbox(['System appears to be STABLE with order q=' num2str(q,numSigDig)], ...
			['Stability test for ''' thisEntry ''''], 'modal');
        else
            msgbox(['System appears to be UNSTABLE with order q=' num2str(q,numSigDig)], ...
			['Stability test for ''' thisEntry ''''], 'modal');
        end
        
        % Set figure name
        set(h, 'NumberTitle', 'off');
        set(h, 'Name', ['''' thisEntry '''' ' ' 'Stability Test Results']);
        
    else
        refreshFotfList(handles);
        errordlg('Object no longer in workspace or invalid!','Error');
    end

% --- Executes on selection change in fotfList.
function fotfList_Callback(hObject, eventdata, handles)
% hObject    handle to fotfList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns fotfList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fotfList


% --- Executes during object creation, after setting all properties.
function fotfList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fotfList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnAdd.
function btnAdd_Callback(hObject, eventdata, handles)
% hObject    handle to btnAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    % Form user prompt
    prompt={'System name:', ...
            'Zero polynomial (e.g. s+1):', ...
            'Pole polynomial (e.g. -s^1.5-1):', ...
            'Delay [sec]:'};
        
    numlines=1;
    
    defaultanswer={'', '', '', '0'};
    name='Create new FO transfer function';
    
    options.Resize      = 'on';
    options.WindowStyle = 'normal';
    
    sysData=inputdlg(prompt,name,numlines,defaultanswer,options);
    
    if ~isempty(sysData)
        sysNewName = sysData{1};

        zeroPoly   = sysData{2};
        polePoly   = sysData{3};
        
        delay      = sysData{4};

        G = newfotf(zeroPoly, polePoly, str2num(delay));

        % Assign in workspace
        assignin('base', sysNewName, G);
        
        % Refresh system list
        refreshFotfList(handles);
        
    end

% --- Executes on button press in btnEdit.
function btnEdit_Callback(hObject, eventdata, handles)
% hObject    handle to btnEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    selected   = get(handles.fotfList, 'Value');
    allEntries = get(handles.fotfList, 'String');
    thisEntry  = allEntries{selected};
    
    [varex, varcl] = varexists(thisEntry);
    
    if varex && varcl
        myVar               = evalin('base', thisEntry);
        [a na b nb ioDel]   = fotfparam(myVar);
        
        % Form user prompt
        prompt={'System name:', ...
                'Zero polynomial (e.g. s+1):', ...
                'Pole polynomial (e.g. -s^1.5-1):', ...
                'Delay [sec]):'};

        numlines=1;
        
		config = fomcon('config');
        numSigDig = config.Core.General.Model_significant_digits;
		
        defaultanswer={thisEntry, fpoly2str(b,nb), fpoly2str(a,na), num2str(ioDel, numSigDig)};
        name='Edit FO transfer function';

        options.Resize      = 'on';
        options.WindowStyle = 'normal';
        
        sysData=inputdlg(prompt,name,numlines,defaultanswer,options);
        
        if ~isempty(sysData)
        
            sysNewName = sysData{1};

            zeroPoly   = sysData{2};
            polePoly   = sysData{3};
            
            delay = sysData{4};

            G = newfotf(zeroPoly, polePoly, str2num(delay));

            % Refresh variable in workspace
            evalin('base', ['clear ' thisEntry]);
            assignin('base', sysNewName, G);
            refreshFotfList(handles);
        
        end

	else
            refreshFotfList(handles);
            errordlg('Object no longer in workspace or invalid!','Error');
	end
    
    
% --- Executes on button press in btnDelete.
function btnDelete_Callback(hObject, eventdata, handles)
% hObject    handle to btnDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Remove list entry
    selected = get(handles.fotfList, 'Value');
    prev_str = get(handles.fotfList, 'String');
    
    thisEntry = false;
    
    if selected
        % Get this entry
        thisEntry = prev_str{selected};
    end
    
    if ~isempty(prev_str)
        prev_str(get(handles.fotfList,'Value')) = [];
        set(handles.fotfList, 'String', prev_str, ...
            'Value', min(selected,length(prev_str)));
    end
    
    if thisEntry
        % Remove workspace variable (if it exists)
        [varex, varcl] = varexists(thisEntry);
        
        if varex && varcl
            evalin('base', ['clear ' thisEntry]);
            refreshFotfList(handles);
        else
            errordlg('Wrong variable class or already deleted!','Error');
        end
    end

% --- Executes on button press in btnBode.
function btnBode_Callback(hObject, eventdata, handles)
% hObject    handle to btnBode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    plotResponse(handles, 'bode');

% --- Executes on button press in btnNyquist.
function btnNyquist_Callback(hObject, eventdata, handles)
% hObject    handle to btnNyquist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    plotResponse(handles, 'nyquist');

% --- Executes on button press in btnGetFreq.
function btnGetFreq_Callback(hObject, eventdata, handles)
% hObject    handle to btnGetFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    fr_resp_types = {'bode', 'nyquist', 'nichols','rootlocus'};
    plotResponse(handles, fr_resp_types{get(handles.menuFreqResponse,'Value')});


function txtFreq_Callback(hObject, eventdata, handles)
% hObject    handle to txtFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFreq as text
%        str2double(get(hObject,'String')) returns contents of txtFreq as a double


% --- Executes during object creation, after setting all properties.
function txtFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnStep.
function btnStep_Callback(hObject, eventdata, handles)
% hObject    handle to btnStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    plotResponse(handles, 'step');

% --- Executes on button press in btnSim.
function btnSim_Callback(hObject, eventdata, handles)
% hObject    handle to btnSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    plotResponse(handles, 'sim');


function txtTime_Callback(hObject, eventdata, handles)
% hObject    handle to txtTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTime as text
%        str2double(get(hObject,'String')) returns contents of txtTime as a double


% --- Executes during object creation, after setting all properties.
function txtTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtU_Callback(hObject, eventdata, handles)
% hObject    handle to txtU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtU as text
%        str2double(get(hObject,'String')) returns contents of txtU as a double


% --- Executes during object creation, after setting all properties.
function txtU_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in btnView.
function btnView_Callback(hObject, eventdata, handles)
% hObject    handle to btnView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    selected = get(handles.fotfList, 'Value');
    allEntries = get(handles.fotfList, 'String');
    thisEntry = allEntries{selected};
    
    [varex, varcl] = varexists(thisEntry);
    
    if varex && varcl
        disp(char(10));
        disp(['FOTF system ''' thisEntry ''':' char(10)]);
        myFotf = evalin('base', thisEntry)
        disp(char(10));
    else
        refreshFotfList(handles);
        errordlg('Object no longer in workspace or invalid!','Error');
    end


% --- Executes on button press in btnRefresh.
function btnRefresh_Callback(hObject, eventdata, handles)
% hObject    handle to btnRefresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    refreshFotfList(handles);
    
function refreshFotfList(handles)
    
    % If systems are present, remember currently selected system
    if ~isempty(get(handles.fotfList,'String'))
        selected   = get(handles.fotfList, 'Value');
        allEntries = get(handles.fotfList, 'String');
        rememberEntry  = allEntries{selected};
    else
        rememberEntry = '';
    end

    % Go through all base workspace variables and populate list
    % if any FOTF systems are present
    allVars  = evalin('base', 'who');
    fotfVars = {};

    for n=1:length(allVars)
       if strcmp(evalin('base',strcat('class(',allVars{n},')')), 'fotf') 
           fotfVars{end+1} = allVars{n};
       end
    end

    if ~isempty(fotfVars)
        % Populate list
        set(handles.fotfList,'String',fotfVars);
        
        % Try to find the remembered system
        toFind = find(strcmp(rememberEntry,fotfVars));        
        if ~isempty(toFind)
            set(handles.fotfList,'Value',toFind);
        else
            set(handles.fotfList,'Value',1);
        end
        
        % Enable buttons
        enableActions(handles, true);
        
    else
        set(handles.fotfList,'String','');
        
        % Disable buttons
        enableActions(handles, false);
        
    end
    
function plotResponse(handles, figureType)
    % Plots desired response type
    selected   = get(handles.fotfList, 'Value');
    allEntries = get(handles.fotfList, 'String');
    thisEntry  = allEntries{selected};
    
    [varex, varcl] = varexists(thisEntry);
    
    if varex && varcl
        myVar   = evalin('base', thisEntry);
        
        % Get values
        % Time vector
		if ~isempty(get(handles.txtTime, 'String'))
			timeVec = evalin('base', get(handles.txtTime, 'String'));
        else
            set(handles.txtTime,'String','0:0.1:30');
            timeVec = 0:0.1:30;
		end
		
		% u(t) vector
		if ~isempty(get(handles.txtU, 'String')) && ...
                length(evalin('base', get(handles.txtU, 'String'))) == length(timeVec)
			uVec = evalin('base', get(handles.txtU, 'String'));      
        else
            set(handles.txtU,'String', ['1*ones(1,' num2str(size(timeVec,2)) ')']);
            uVec = 1*ones(1,size(timeVec,2));
        end
        
        % Frequency range
        if ~isempty(get(handles.txtFreq, 'String'))
            freqRange = evalin('base', get(handles.txtFreq, 'String'));
        else
            set(handles.txtFreq,'String','logspace(-4,4,5000)');
            freqRange = logspace(-4,4,5000);
        end
        
        % Retrieve stored GUI data
        gui_data = get(handles.output, 'UserData');
        time_h   = gui_data.Fig_time_response;
        freq_h   = gui_data.Fig_freq_response;
        
        % Check whether to superimpose time domain figures
        if get(handles.chkTimeSuperimpose, 'Value')
            if isempty(time_h)
                time_h = figure;
                gui_data.Fig_time_response = time_h;
            else
                time_h = figure(time_h);
            end
            h = time_h;
            hold on;
        elseif get(handles.chkFreqSuperimpose, 'Value')
            if isempty(freq_h)
                freq_h = figure;
                gui_data.Fig_freq_response = freq_h;
            else
                freq_h = figure(freq_h);
            end
            h = freq_h;
            hold on;
        else
            h  = figure;
        end

        % Save figure data
        set(handles.output, 'UserData', gui_data);
        
        yy = [];
        
        % Plot requested figure
        switch figureType
            case 'step'
                figureName = 'Step Response';
                yy = step(myVar, timeVec);
                plot(timeVec, yy);
            case 'impulse'
                figureName = 'Impulse Response';
                yy = impulse(myVar, timeVec);
                plot(timeVec, yy);
            case 'sim'
                figureName = 'Time Response Simulation';
                y=lsim(myVar,uVec, timeVec);
                plot(timeVec, y);
                yy = y;
            case 'bode'
                figureName = 'Bode Diagram';
                bode(myVar, freqRange);
            case 'nyquist'
                figureName = 'Nyquist Diagram';
                nyquist(myVar, freqRange);
            case 'nichols'
                figureName = 'Nichols Diagram';
                nichols(myVar, freqRange);
            case 'rootlocus'
                figureName = 'Root Locus';
                rlocus(myVar);
            otherwise
                error('Unknown plot type requested');
        end
        
        % Save data, if requested
        ySaveName = get(handles.txtSaveY, 'String');
        if ~isempty(yy) && ~isempty(ySaveName)
            assignin('base', ySaveName, yy)
        end
        
        % Set axes properties for time-domain simulation
        switch figureType
            case {'step', 'impulse', 'sim'}
                xlabel('Time [sec]');
                ylabel('Amplitude');
        end
				
        % Set figure name
        set(h, 'NumberTitle', 'off');
        set(h, 'Name', ['''' thisEntry '''' ' ' figureName]);
        
        % Turn grid on
        grid;
        
        % Go through axes
        ax = get(h, 'Children');
        for n=1:length(ax)
            
            if strcmp(get(ax(n), 'Type'), 'axes')
                % Set grids
                set(ax(n), 'XGrid', 'on', 'YGrid', 'on');
            end
            
        end
        
    else
        refreshFotfList(handles.fotfList);
        errordlg('Object no longer in workspace or invalid!','Error');
    end


% --------------------------------------------------------------------
function menuTools_Callback(hObject, eventdata, handles)
% hObject    handle to menuTools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuFracPid_Callback(hObject, eventdata, handles)
% hObject    handle to menuFracPid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    selected   = get(handles.fotfList, 'Value');
    allEntries = get(handles.fotfList, 'String');
    thisEntry  = allEntries{selected};
    
    [varex, varcl] = varexists(thisEntry);
    
    if varex && varcl
        fpid('UserData', thisEntry);
    else
        refreshFotfList(handles.fotfList);
        errordlg('Object no longer in workspace or invalid!','Error');
    end


% --- Executes on selection change in listExport.
function listExport_Callback(hObject, eventdata, handles)
% hObject    handle to listExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listExport contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listExport

% Enable/disable 'launch LTI viewer option'
switch(get(hObject,'Value'))
    case {1, 2}
        set(handles.chkLaunchLTIView, 'Enable', 'on');
    otherwise
        set(handles.chkLaunchLTIView, 'Enable', 'off');
end

% --- Executes during object creation, after setting all properties.
function listExport_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnExport.
function btnExport_Callback(hObject, eventdata, handles)
% hObject    handle to btnExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    % Check FOTF
    selected   = get(handles.fotfList, 'Value');
    allEntries = get(handles.fotfList, 'String');
    thisEntry  = allEntries{selected};
    
    [varex, varcl] = varexists(thisEntry);
    
    if ~(varex && varcl)
        refreshFotfList(handles.fotfList);
        error('Object no longer in workspace or invalid!');
    end
    
    % Get FOTF
    myG = evalin('base', thisEntry);

    exportType = get(handles.listExport, 'Value');
    
    % Export data based on pop-up menu
    switch(exportType)
        
        case {1, 2}
            % Oustaloup filter zpk
            prompt={'ZPK workspace name:', ...
            'Low frequency bound wb [rad/s]', ...
            'High frequency bound wh [rad/s]', ...
            'Order of approximation:'};
            
            if exportType == 1
                name='Export as Oustaloup Filter ZPK';
            else
                name='Export as Refined Oustaloup Filter ZPK';
            end
        
            numlines=1;

            defaultAnswer = {'', '0.0001', '10000', '5'};

            options.Resize      = 'on';
            options.WindowStyle = 'normal';
            
            ofData=inputdlg(prompt,name,numlines,defaultAnswer,options);
            
            if ~isempty(ofData)
            
                ofName  = ofData{1};
                ofWb    = str2num(ofData{2});
                ofWh    = str2num(ofData{3});
                ofN     = str2num(ofData{4});

                if exportType == 1
                    
                    try
                        exportModel = oustapp(myG,ofWb,ofWh,ofN,'oust');
                        assignin('base', ofName, exportModel);
                    catch
                        % Bad model
                        errordlg('Erroneous model, export aborted: possibly due to too high order of approximation', 'Export failed');
                    end
                    
                elseif exportType == 2
                    
                    try
                        exportModel = oustapp(myG,ofWb,ofWh,ofN,'ref');
                        assignin('base', ofName, exportModel);
                    catch
                        errordlg('Erroneous model, export aborted: possibly due to too high order of approximation', 'Export failed');
                    end
                    
                end
                
                % Launch LTI viewer if required
                if (get(handles.chkLaunchLTIView,'Value'))
                   evalin('base', ['ltiview(' ofName ')']);
                end
            
            end
            
            
        case 3
            % Fractional State-space
            options.Resize      = 'on';
            options.WindowStyle = 'normal';
            
            ssData = inputdlg({'Fractional state-space workspace name:'},...
                              'Export fractional state-space', ...
                              1, ...
                              {''}, ...
                              options);
            
            if ~isempty(ssData)                     
                assignin('base', ssData{1}, tf2ss(myG));
            end
        
        case 4
            
            % CRONE frac_tf
            options.Resize      = 'on';
            options.WindowStyle = 'normal';
            
            ssData = inputdlg({'CRONE frac_tf workspace name:'},...
                              'Export CRONE frac_tf', ...
                              1, ...
                              {''}, ...
                              options);
            
            if ~isempty(ssData)
                assignin('base', ssData{1}, tf2tf_c(myG));
            end
            
        case 5
            
            % CRONE frac_ss
            options.Resize      = 'on';
            options.WindowStyle = 'normal';
            
            ssData = inputdlg({'CRONE frac_ss workspace name:'},...
                              'Export CRONE frac_tf', ...
                              1, ...
                              {''}, ...
                              options);
            
            if ~isempty(ssData)
                assignin('base', ssData{1}, tf2ss_c(myG));
            end
            
        otherwise
            
            errordlg('Action undefined!','Error');
    end

% --- Executes on button press in btnInputFreq.
function btnInputFreq_Callback(hObject, eventdata, handles)
% hObject    handle to btnInputFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

            prompt={'Lower frequency [rad/s]:', ...
            'Higher frequency [rad/s]:', ...
            'Data points:'};
            
            name='Input frequency parameters';
        
            numlines=1;

            defaultAnswer = {'0.001','1000','5000'};

            options.Resize      = 'on';
            options.WindowStyle = 'normal';
            
            ftData=inputdlg(prompt,name,numlines,defaultAnswer,options);
            
            if ~isempty(ftData)
            
                lFreq = get10exp(str2num(ftData{1}));
                hFreq = get10exp(str2num(ftData{2}));

                freq = ['logspace(' num2str(lFreq) ',' num2str(hFreq) ',' ftData{3} ')'];
                set(handles.txtFreq, 'String', freq);
            
            end


% --- Executes on button press in chkShowProgBar.
function chkShowProgBar_Callback(hObject, eventdata, handles)
% hObject    handle to chkShowProgBar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of chkShowProgBar
    
    % Enable/disable progress bar
    switch get(hObject, 'Value')
        case 0
            show_pbar off
        case 1
            show_pbar on
    end
        

% --- Executes on button press in btnTimeSim.
function btnTimeSim_Callback(hObject, eventdata, handles)
% hObject    handle to btnTimeSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resp_types = {'step', 'impulse', 'sim'};
plotResponse(handles, resp_types{get(handles.menuTimeResponse, 'Value')});

% --- Executes on button press in btnUtInput.
function btnUtInput_Callback(hObject, eventdata, handles)
% hObject    handle to btnUtInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
            
            prompt={'Start time [0, ...):', ...
            'End time:', ...
            'Time step:', ...
            'Set value:'};
            
            name='Input simulation parameters';
        
            numlines=1;

            defaultAnswer = {'0', '30', '0.1', '1'};

            options.Resize      = 'on';
            options.WindowStyle = 'normal';
            
            utData=inputdlg(prompt,name,numlines,defaultAnswer,options);
            
            if ~isempty(utData)
            
                timeVecS = [utData{1} ':' utData{3} ':' utData{2}];
                set(handles.txtTime, 'String', timeVecS);

                timeVec = evalin('base', get(handles.txtTime, 'String'));
                set(handles.txtU, 'String', [utData{4} '*ones(1,' num2str(size(timeVec,2)) ')']);
            
            end

% --- Executes on button press in btnSaveUt.
function btnSaveUt_Callback(hObject, eventdata, handles)
% hObject    handle to btnSaveUt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
            
            prompt={'u(t) workspace name', ...
                    't workspace name:', ...
                   };
            
            name='Export simulation initial data';
        
            numlines=1;

            defaultAnswer = {'u', 't'};

            options.Resize      = 'on';
            options.WindowStyle = 'normal';
            
            utData=inputdlg(prompt,name,numlines,defaultAnswer,options);
            
            if ~isempty(utData)
                assignin('base', utData{1}, str2num(get(handles.txtU, 'String'))');
                assignin('base', utData{2}, str2num(get(handles.txtTime, 'String'))');
            end

% --- Enables all buttons/actions
function enableActions(handles, status)

    if status
        enableElement = 'on';
    else
        enableElement = 'off';
    end
    
    % Enable/disable UI elements
    set(handles.btnEdit,'Enable',enableElement);    % Edit
    set(handles.btnDelete,'Enable',enableElement);  % Delete
    
    set(handles.btnExport,'Enable',enableElement);  % Export
    
    set(handles.btnView,'Enable',enableElement);    % View
    set(handles.btnStable,'Enable',enableElement);  % Stability test
    
    set(handles.menuTimeResponse,'Enable',enableElement);
    set(handles.btnTimeSim,'Enable',enableElement);     % Simulate
    
    set(handles.menuFreqResponse,'Enable',enableElement);
    set(handles.btnGetFreq,'Enable',enableElement);    % Bode
    
    set(handles.listExport,'Enable',enableElement);    % Export format list
    
    set(handles.menuFracPid,'Enable',enableElement); % Menu:Fractional PID Tool


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuTDIdent_Callback(hObject, eventdata, handles)
% hObject    handle to menuTDIdent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in chkLaunchLTIView.
function chkLaunchLTIView_Callback(hObject, eventdata, handles)
% hObject    handle to chkLaunchLTIView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkLaunchLTIView


% --------------------------------------------------------------------
function menuTimeDomainId_Callback(hObject, eventdata, handles)
% hObject    handle to menuTimeDomainId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Open FOTFID tool
    fotfid();

% --------------------------------------------------------------------
function menuFreqDomainId_Callback(hObject, eventdata, handles)
% hObject    handle to menuFreqDomainId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    % Open FOTFRID tool
    fotfrid();


% --- Executes on selection change in menuTimeResponse.
function menuTimeResponse_Callback(hObject, eventdata, handles)
% hObject    handle to menuTimeResponse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menuTimeResponse contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuTimeResponse


% --- Executes during object creation, after setting all properties.
function menuTimeResponse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuTimeResponse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menuFreqResponse.
function menuFreqResponse_Callback(hObject, eventdata, handles)
% hObject    handle to menuFreqResponse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menuFreqResponse contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuFreqResponse


% --- Executes during object creation, after setting all properties.
function menuFreqResponse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuFreqResponse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtSaveY_Callback(hObject, eventdata, handles)
% hObject    handle to txtSaveY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtSaveY as text
%        str2double(get(hObject,'String')) returns contents of txtSaveY as a double


% --- Executes during object creation, after setting all properties.
function txtSaveY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSaveY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menuConfig_Callback(hObject, eventdata, handles)
% hObject    handle to menuConfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fomcon('config');


% --- Executes on button press in chkFreqSuperimpose.
function chkFreqSuperimpose_Callback(hObject, eventdata, handles)
% hObject    handle to chkFreqSuperimpose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkFreqSuperimpose
if get(handles.chkFreqSuperimpose, 'Value')
    set(handles.chkTimeSuperimpose, 'Value', 0);
end


% --- Executes on button press in chkTimeSuperimpose.
function chkTimeSuperimpose_Callback(hObject, eventdata, handles)
% hObject    handle to chkTimeSuperimpose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkTimeSuperimpose
if get(handles.chkTimeSuperimpose, 'Value')
    set(handles.chkFreqSuperimpose, 'Value', 0);
end
