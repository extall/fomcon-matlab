function varargout = fpid(varargin)
% FPID Fractional-order PID controller design tool.
%
% This tool implements guided fractional-order PID controller design
% through integer-order methods and optimization.
%
% Last Modified by GUIDE v2.5 09-Oct-2013 17:13:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fpid_OpeningFcn, ...
                   'gui_OutputFcn',  @fpid_OutputFcn, ...
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


% --- Executes just before fpid is made visible.
function fpid_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fpid (see VARARGIN)

% Get and set initial data (system name)
sysName  = get(hObject, 'UserData');
if ~isempty(sysName)
    set(handles.txtSystem, 'String', sysName);
end

% Set time and step values
set(handles.txtTime, 'String', '0:0.1:100');
set(handles.txtSV, 'String', '1');

% Set PID parameters
set(handles.txtKp, 'String', '1');
set(handles.txtKi, 'String', '1');
set(handles.txtLambda, 'String', '0.5');
set(handles.txtKd, 'String', '1');
set(handles.txtMu, 'String', '0.5');

% Choose default command line output for fpid
handles.output = hObject;

% Put image inside of picture axes
pidPlant = importdata('fpid.jpg');
axes(handles.axSchematic);
hi = imagesc(pidPlant);
axis off;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fpid wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = fpid_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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


% --- Executes on button press in btnView.
function btnView_Callback(hObject, eventdata, handles)
% hObject    handle to btnView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    sysName = get(handles.txtSystem, 'String');
    [varex, varcl, varclass] = varexists(sysName);
    
    if varex && varcl
        % View fractional order controller and system information
        
        % Fetch PID parameters
        [Kp,Ki,lambda,Kd,mu] = getPidParams(handles);
    
        % Get corrsponding FOTF PID controller and its type
        [myPid, pidType] = fracpid(Kp,Ki,lambda,Kd,mu);
        
        % Get plant
        myPlant = evalin('base', sysName);
        
        % Full control system
        cPair    = myPid * myPlant;
        fullCtrl = feedback(cPair, 1);
        
        % Controller information
        disp(char(13));
        disp(['Current Controller (type: fractional ' pidType '):']);
        myPid
        
        % System information
        disp(char(13));
        disp('Current Plant:');
        myPlant
        
        % Full control system
        disp(char(13));
        disp('Full control system:');
        fullCtrl
        disp(char(13));
    elseif varex
        
        switch(varclass)
            case {'tf', 'zpk', 'ss'}
                % Fetch PID parameters
                [Kp,Ki,lambda,Kd,delta] = getPidParams(handles);
    
                % Get corrsponding FOTF PID controller and its type
                [myPid, pidType] = fracpid(Kp,Ti,lambda,Kd,mu);
        
                % Controller information
                disp(['Current Controller (type: ' pidType '):']);
                myPid
                
            otherwise
                errordlg('Object is not a LTI model!', 'Error');
        end
        
    else
        set(handles.txtSystem, 'String', '');
        errordlg('Object no longer in workspace or invalid!', 'Error');
    end    

% --- Executes on button press in btnSimulate.
function btnSimulate_Callback(hObject, eventdata, handles)
% hObject    handle to btnSimulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    sysName = get(handles.txtSystem, 'String');
    [varex, varcl, varclass] = varexists(sysName);
    
    % Time vector
	if ~isempty(get(handles.txtTime, 'String'))
		timeVec = evalin('base', get(handles.txtTime, 'String'));
    else
        set(handles.txtTime, 'String', '0:0.1:30');
        timeVec = 0:0.1:30;
	end
		
	% SV
	if ~isempty(get(handles.txtSV, 'String'))
        SVval = evalin('base', get(handles.txtSV, 'String')); 
        SV = SVval*ones(1,size(timeVec,2));
	else
        set(handles.txtSV, 'String', '1');
        SV = ones(1,size(timeVec,2));
    end
    
    % Check for delay system
    if varex
        plant = evalin('base', sysName);
        if ~fleq(plant.ioDelay, 0)
            hasDelay = 1;
        else
            hasDelay = 0;
        end
    end
    
    if varex && varcl && ~hasDelay
        
        % Fetch PID parameters
        [Kp,Ki,lambda,Kd,mu] = getPidParams(handles);
    
        % Get corrsponding FOTF PID controller and its type
        [myPid, pidType] = fracpid(Kp,Ki,lambda,Kd,mu);
        
        % Get plant
        myPlant = evalin('base', sysName);
        
        % Full control system
        cPair    = myPid * myPlant;
        fullCtrl = feedback(cPair,1);
        
        % Setup figure
        h  = figure;
        
        % Plot figure
        figureName = [pidType ' Control System Time Response Simulation'];
        y=lsim(fullCtrl, SV, timeVec);
        plot(timeVec, y);
        
        % Plot SV line
        line([timeVec(1) timeVec(size(timeVec,2))], ...
                [SV(1) SV(1)], 'Color', 'red',  ...
                'LineWidth', 1, 'LineStyle', '-');
        
        
        % Set figure name
        set(h, 'NumberTitle', 'off');
        set(h, 'Name', ['''' sysName '''' ' ' figureName]);
        
        % Go through axes
        ax = get(h, 'Children');
        for n=1:length(ax)
            
            if strcmp(get(ax(n), 'Type'), 'axes')
                % Set grids
                set(ax(n), 'XGrid', 'on', 'YGrid', 'on');
            end
            
        end
        
    elseif varex || (varcl && hasDelay)
        
        switch(varclass)
            
            case {'tf', 'zpk', 'ss', 'fotf'}
                
                % Simulation parameters
                name='PID approximation parameters';
                
                prompt={'Approximation type (oust or ref):', ...
                'Low frequency bound wb [rad/s]', ...
                'High frequncy bound wh [rad/s]', ...
                'Order of approximation:'};

                numlines=1;

                defaultAnswer = {'oust', '0.0001', '10000', '5'};

                options.Resize      = 'on';
                options.WindowStyle = 'normal';

                ofData=inputdlg(prompt,name,numlines,defaultAnswer,options);
            
                if ~isempty(ofData)
            
                    ofType  = ofData{1};
                    ofWb    = str2num(ofData{2});
                    ofWh    = str2num(ofData{3});
                    ofN     = str2num(ofData{4});
                    
                    % Get plant
                    myPlant = evalin('base', sysName);
                    
                    % Fetch PID parameters
                    [Kp,Ki,lambda,Kd,mu] = getPidParams(handles);
    
                    % Get corrsponding FOTF PID controller and its type
                    [myPid, pidType] = fracpid(Kp,Ki,lambda,Kd,mu);
                    
                    % PID approximation
                    pidApprox = oustapp(myPid, ofWb, ofWh, ofN, ofType);
                    
                    if ~isproper(pidApprox)
                        pidApprox = toproper(pidApprox, ofWh);
                    end
                    
                    % If plant is given by fotf, convert with same params
                    if strcmp(varclass, 'fotf')
                        myPlant = oustapp(myPlant, ofWb, ofWh, ofN, ofType);
                    end
                    
                    % Convert plant
                    myPlant = ss(myPlant);
                    
                    % Get control system
                    ctrlSys = feedback(pidApprox*myPlant, 1);
                    assignin('base','gg',ctrlSys);
                    
                    % Setup figure
                    h  = figure;

                    % Plot figure
                    figureName = [pidType ' Control System Time Response Simulation'];
                    y=lsim(ctrlSys, SV, timeVec);
                    plot(timeVec, y);

                    % Plot SV line
                    line([timeVec(1) timeVec(size(timeVec,2))], ...
                            [SV(1) SV(1)], 'Color', 'red',  ...
                            'LineWidth', 1, 'LineStyle', '-');


                    % Set figure name
                    set(h, 'NumberTitle', 'off');
                    set(h, 'Name', ['''' sysName '''' ' ' figureName]);

                    % Go through axes
                    ax = get(h, 'Children');
                    for n=1:length(ax)

                        if strcmp(get(ax(n), 'Type'), 'axes')
                            % Set grids
                            set(ax(n), 'XGrid', 'on', 'YGrid', 'on');
                        end

                    end
                
                end
                
            otherwise
                errordlg('Object is not a valid LTI model!', 'Error');
        end
        
    else
        set(handles.txtSystem, 'String', '');
        errordlg('Object no longer in workspace or invalid!', 'Error');
    end    

function txtKp_Callback(hObject, eventdata, handles)
% hObject    handle to txtKp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtKp as text
%        str2double(get(hObject,'String')) returns contents of txtKp as a double


% --- Executes during object creation, after setting all properties.
function txtKp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtKp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtKi_Callback(hObject, eventdata, handles)
% hObject    handle to txtKi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtKi as text
%        str2double(get(hObject,'String')) returns contents of txtKi as a double


% --- Executes during object creation, after setting all properties.
function txtKi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtKi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtLambda_Callback(hObject, eventdata, handles)
% hObject    handle to txtLambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtLambda as text
%        str2double(get(hObject,'String')) returns contents of txtLambda as a double


% --- Executes during object creation, after setting all properties.
function txtLambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtLambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function txtKd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtKd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtMu_Callback(hObject, eventdata, handles)
% hObject    handle to txtMu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtMu as text
%        str2double(get(hObject,'String')) returns contents of txtMu as a double


% --- Executes during object creation, after setting all properties.
function txtMu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtMu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtSystem_Callback(hObject, eventdata, handles)
% hObject    handle to txtSystem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtSystem as text
%        str2double(get(hObject,'String')) returns contents of txtSystem as a double


% --- Executes during object creation, after setting all properties.
function txtSystem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSystem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtSV_Callback(hObject, eventdata, handles)
% hObject    handle to txtSV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtSV as text
%        str2double(get(hObject,'String')) returns contents of txtSV as a double


% --- Executes during object creation, after setting all properties.
function txtSV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end        

% --- Executes on button press in btnExportPID.
function btnExportPID_Callback(hObject, eventdata, handles)
% hObject    handle to btnExportPID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    % Prompt user for workspace variable to export to
    name='Export PID-FOTF to Workspace';
    prompt={'Workspace variable name:'};
    defaultanswer = {''};
    numlines=1;
    options.WindowStyle='normal';
    toExport=inputdlg(prompt,name,numlines,defaultanswer,options);
    
    if ~isempty(toExport)
        
        varName =toExport{1};

        % Fetch PID parameters
        [Kp,Ki,lambda,Kd,mu] = getPidParams(handles);

        % Get corrsponding FOTF PID controller
        myPid = fracpid(Kp,Ki,lambda,Kd,mu);

        % Export PID controller FOTF object
        assignin('base', varName, myPid);
        
    end
    
% --- Fetches PID parameters
function [Kp,Ki,lambda,Kd,mu] = getPidParams(handles)
    
    % Get all values from textboxes
    tKp     = get(handles.txtKp, 'String');
    tKi     = get(handles.txtKi, 'String');
    tLambda = get(handles.txtLambda, 'String');
    tKd     = get(handles.txtKd, 'String');
    tMu  = get(handles.txtMu, 'String');
    
    % Check if any of the parameters are empty, use default values for
    % those that are indeed empty
    if isempty(tKp)
        set(handles.txtKp, 'String', '1');
        Kp = 1;
    else
        Kp = evalin('base', tKp);
    end
    
    if isempty(tKi)
        set(handles.txtKi, 'String', '1');
        Ki = 1;
    else
        Ki = evalin('base', tKi);
    end
    
    if isempty(tLambda)
        set(handles.txtLambda, 'String', '0.5');
        lambda = 1;
    else
        lambda = evalin('base', tLambda);
    end
    
    if isempty(tKd)
        set(handles.txtKd, 'String', '0');
        Kd = 0;
    else
        Kd = evalin('base', tKd);
    end
    
    if isempty(tMu)
        set(handles.txtMu, 'String', '0.5');
        mu = 0.5;
    else
        mu = evalin('base', tMu);
    end
        


% --- Executes during object creation, after setting all properties.
function axSchematic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axSchematic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axSchematic


% --------------------------------------------------------------------
function menuTuning_Callback(hObject, eventdata, handles)
% hObject    handle to menuTuning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuIopid_Callback(hObject, eventdata, handles)
% hObject    handle to menuIopid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
% function menuZxc_Callback(hObject, eventdata, handles)
% % hObject    handle to menuZxc (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%     
%     % Get data
%     fotfPlant = get(handles.txtSystem, 'String');
%     [varex, varcl] = varexists(fotfPlant);
%     
%     if varex && varcl
%         
%         % Get plant model
%         fotfPlant = evalin('base', fotfPlant);
%         
%         % Check object
%         [a,na,b,nb] = fotfparam(fotfPlant);
%         
%         numer = (length(b) == 1) && (b(1) == 1) && (nb(1) == 0); 
%         denom = (length(a) == 3) && (na(3) == 0);
%     
%         % Restrict object class
%         if ~(numer && denom)
%             errordlg('Cannot use Zhao-Xue-Chen tuning method for this class of object', 'Error');
%             return;
%         end
%         
%         % Open the tuning utility
%         zxc_gui('UserData', get(handles.txtSystem,'String'));
%         
%     else
%        
%         errordlg('Object no longer in workspace or invalid!', 'Error');
%         
%     end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in btnExportControl.
function btnExportControl_Callback(hObject, eventdata, handles)
% hObject    handle to btnExportControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    sysName = get(handles.txtSystem, 'String');
    [varex, varcl, varclass] = varexists(sysName);
    
    % Fetch PID parameters
    [Kp,Ki,lambda,Kd,mu] = getPidParams(handles);
    
    % Get corrsponding FOTF PID controller and its type
    [myPid, pidType] = fracpid(Kp,Ki,lambda,Kd,mu);
    
    if varex && varcl
        
        % Get plant
        myPlant = evalin('base', sysName);
        
        % Full control system
        cPair    = myPid * myPlant;
        fullCtrl = feedback(cPair,1);
        
        % Prompt user for workspace variable to export to
        name='Save control system';
        prompt={'Workspace variable name:'};
        defaultanswer = {[sysName '_control']};
        numlines=1;
        options.WindowStyle='normal';
        options.Resize = 'on';
        toExport=inputdlg(prompt,name,numlines,defaultanswer,options);
        
        if ~isempty(toExport)
            
            % Get export variable name
            ctrlSys = toExport{1};

            % Save control system
            if ~strcmp(ctrlSys, '')
                assignin('base', ctrlSys, fullCtrl);
            end
            
        end
        
    elseif varex
        
        switch(varclass)
            
            case {'tf', 'zpk', 'ss'}
                
                % Parameters
                name='Export parameters';
                
                prompt={'Workspace variable name:', ...
                'Approximation type (oust or ref):', ...
                'Low frequency bound wb [rad/s]', ...
                'High frequncy bound wh [rad/s]', ...
                'Order of approximation:'};

                numlines=1;

                defaultAnswer = {[sysName '_control'], 'oust', '0.0001', '10000', '5'};

                options.Resize      = 'on';
                options.WindowStyle = 'normal';

                ofData=inputdlg(prompt,name,numlines,defaultAnswer,options);
            
                if ~isempty(ofData)
                    
                    ctrlSys = ofData{1};
                    ofType  = ofData{2};
                    ofWb    = str2num(ofData{3});
                    ofWh    = str2num(ofData{4});
                    ofN     = str2num(ofData{5});
                    
                    % Get plant
                    myPlant = evalin('base', sysName);
                    
                    % Fetch PID parameters
                    [Kp,Ki,lambda,Kd,mu] = getPidParams(handles);
    
                    % Get corrsponding FOTF PID controller and its type
                    [myPid, pidType] = fracpid(Kp,Ki,lambda,Kd,mu);
                    
                    % PID approximation
                    pidApprox = oustapp(myPid, ofWb, ofWh, ofN, ofType);
                    
                    % Convert plant
                    myPlant = ss(myPlant);
                    
                    % Get control system
                    fullCtrl = feedback(pidApprox*myPlant, 1);
                   
                    assignin('base', ctrlSys, fullCtrl);
                    
                end
                
            otherwise
                errordlg('Object is not a valid LTI model!', 'Error');
        end

    else
        set(handles.txtSystem, 'String', '');
        errordlg('Object no longer in workspace or invalid!', 'Error');
    end


% --------------------------------------------------------------------
function menuOptimize_Callback(hObject, eventdata, handles)
% hObject    handle to menuOptimize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    % Get data
    fotfPlant = get(handles.txtSystem, 'String');
    [varex, varcl, varclass] = varexists(fotfPlant);
    
    if varex && (varcl || strcmpi(varclass, 'tf') || ...
                          strcmpi(varclass, 'zpk') || ...
                          strcmpi(varclass, 'ss'))
        
        % Generate optimization utility object
        toOpt.G = get(handles.txtSystem,'String');
        toOpt.Kp = get(handles.txtKp,'String');
        toOpt.Ki = get(handles.txtKi,'String');
        toOpt.Kd = get(handles.txtKd,'String');
        toOpt.lam = get(handles.txtLambda,'String');
        toOpt.mu = get(handles.txtMu,'String');
        
        % Open the optimization utility
        fpid_optim('UserData', toOpt);
        
    else
        
        errordlg('Object no longer in workspace or invalid!', 'Error');
        
    end


% --------------------------------------------------------------------
function menuIntegerPID_Callback(hObject, eventdata, handles)
% hObject    handle to menuIntegerPID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Get data
    fotfPlant = get(handles.txtSystem, 'String');
    [varex, varcl] = varexists(fotfPlant);
    
    if varex && varcl
        
        % Open the integer-order tuning tool
        iopid_tune('UserData', get(handles.txtSystem,'String'));
        
	else
        
        errordlg('Object no longer in workspace or invalid!', 'Error');
        
    end


% --------------------------------------------------------------------
function menuImport_Callback(hObject, eventdata, handles)
% hObject    handle to menuImport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuImportPID_Callback(hObject, eventdata, handles)
% hObject    handle to menuImportPID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    pr = inputdlg({'PID workspace name:'}, ...
                   'Import fractional or integer-order PID', ...
                   1, ...
                   {''});

    if ~isempty(pr)
       modelName = pr{1};
       myModel = evalin('base', modelName);
       [Kp, Ki, ilam, Kd, dmu] = getpid(myModel);
       
	   config = fomcon('config');
       numSigDig = config.Core.General.Model_significant_digits;
	   
       % Set gains/exponents
       set(handles.txtKp,     'String', num2str(Kp,numSigDig));
       set(handles.txtKi,     'String', num2str(Ki,numSigDig));
       set(handles.txtLambda, 'String', num2str(ilam,numSigDig));
       set(handles.txtKd,     'String', num2str(Kd,numSigDig));
       set(handles.txtMu,  'String', num2str(dmu,numSigDig));
       
    end


% --- Executes on button press in btnRealize.
function btnRealize_Callback(hObject, eventdata, handles)
% hObject    handle to btnRealize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        toImp.Kp  = get(handles.txtKp,'String');
        toImp.Ki  = get(handles.txtKi,'String');
        toImp.Kd  = get(handles.txtKd,'String');
        toImp.lam = get(handles.txtLambda,'String');
        toImp.mu  = get(handles.txtMu,'String');
        
        impid('UserData', toImp);


% --- Executes on button press in btnOLBode.
function btnOLBode_Callback(hObject, eventdata, handles)
% hObject    handle to btnOLBode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    sysName = get(handles.txtSystem, 'String');
    [varex, varcl, varclass] = varexists(sysName);
    
    % Fetch PID parameters
    [Kp,Ki,lambda,Kd,mu] = getPidParams(handles);
    
    % Get corrsponding FOTF PID controller and its type
    [myPid, pidType] = fracpid(Kp,Ki,lambda,Kd,mu);

    if varex && varcl
        
        % Get plant
        myPlant = evalin('base', sysName);
        
        % Full control system
        cPair    = myPid * myPlant;
        
        % Parameters
        name='Bode plot parameters';
        
        prompt={'Frequencies of interest [rad/s]'};
        
        numlines=1;
        
        defaultAnswer = {'logspace(-5,5,1000)'};
        
        options.Resize      = 'on';
        options.WindowStyle = 'normal';
        
        ofData=inputdlg(prompt,name,numlines,defaultAnswer,options);
        
        if ~isempty(ofData)
            
            w_ev    = evalin('base', ofData{1});
            
            h = figure;
            bode(cPair, w_ev);
            grid;
            set(h, 'Name', [pidType ' * ' sysName ' open-loop frequency domain response']);
        end
        
    elseif varex
        
        switch(varclass)
            
            case {'tf', 'zpk', 'ss'}
                
                % Parameters
                name='Bode plot parameters';
                
                prompt={'Approximation type (oust or ref):', ...
                'Low frequency bound wb [rad/s]', ...
                'High frequncy bound wh [rad/s]', ...
                'Order of approximation:', ...
                'Frequencies of interest [rad/s]'};

                numlines=1;

                defaultAnswer = {'oust', '0.0001', '10000', '5', 'logspace(-5,5,1000)'};

                options.Resize      = 'on';
                options.WindowStyle = 'normal';

                ofData=inputdlg(prompt,name,numlines,defaultAnswer,options);
            
                if ~isempty(ofData)
                    
                    ofType  = ofData{1};
                    ofWb    = str2num(ofData{2});
                    ofWh    = str2num(ofData{3});
                    ofN     = str2num(ofData{4});
                    w_ev    = evalin('base', ofData{5});
                    
                    % Get plant
                    myPlant = evalin('base', sysName);
                    
                    % Fetch PID parameters
                    [Kp,Ki,lambda,Kd,mu] = getPidParams(handles);
    
                    % Get corrsponding FOTF PID controller and its type
                    [myPid, pidType] = fracpid(Kp,Ki,lambda,Kd,mu);
                    
                    % PID approximation
                    pidApprox = oustapp(myPid, ofWb, ofWh, ofN, ofType);
                    
                    % Convert plant
                    myPlant = ss(myPlant);
                    
                    % Get control system
                    cPair = pidApprox*myPlant;
                    
                    h = figure;
                    bode(cPair, w_ev);
                    grid;
                    set(h, 'Name', [pidType ' * ' sysName ' open-loop frequency domain response']);
                end
                
            otherwise
                errordlg('Object is not a valid LTI model!', 'Error');
        end

    else
        set(handles.txtSystem, 'String', '');
        errordlg('Object no longer in workspace or invalid!', 'Error');
    end


% --------------------------------------------------------------------
function menuImportPIDConfig_Callback(hObject, eventdata, handles)
% hObject    handle to menuImportPIDConfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,path] = uigetfile();
if ~isequal(filename, 0) && ~isequal(path, 0)
    config_struct   = load([path filename]);
    config = config_struct.FPID_Optimizer_GUI_config;
    
    % Handles shortcut
    h = handles;
    FPIDParams = config.FPIDParams;
    set(h.txtKp, 'String', FPIDParams.Kp);
    set(h.txtKi, 'String', FPIDParams.Ki);
    set(h.txtKd, 'String', FPIDParams.Kd);
    set(h.txtLambda, 'String', FPIDParams.Lam);
    set(h.txtMu, 'String', FPIDParams.Mu);
end



function txtKd_Callback(hObject, eventdata, handles)
% hObject    handle to txtKd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtKd as text
%        str2double(get(hObject,'String')) returns contents of txtKd as a double
