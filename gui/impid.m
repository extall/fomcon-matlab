function varargout = impid(varargin)
% IMPID Implement a PID controller
%
% This tool is used to implement digital (fractional-order) PID controllers
% based on frequency-domain analysis.
%
% Last Modified by GUIDE v2.5 20-Jun-2012 02:48:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @impid_OpeningFcn, ...
                   'gui_OutputFcn',  @impid_OutputFcn, ...
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


% --- Executes just before impid is made visible.
function impid_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to impid (see VARARGIN)

% Get implementation data, if supplied
toImp  = get(hObject, 'UserData');

if ~isempty(toImp)
    set(handles.txtKp, 'String', toImp.Kp);
    set(handles.txtKi, 'String', toImp.Ki);
    set(handles.txtKd, 'String', toImp.Kd);
    set(handles.txtLam, 'String', toImp.lam);
    set(handles.txtMu, 'String', toImp.mu);
end

% Choose default command line output for impid
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes impid wait for user response (see UIRESUME)
% uiwait(handles.impid_gui);


% --- Outputs from this function are returned to the command line.
function varargout = impid_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function txtDName_Callback(hObject, eventdata, handles)
% hObject    handle to txtDName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtDName as text
%        str2double(get(hObject,'String')) returns contents of txtDName as a double


% --- Executes during object creation, after setting all properties.
function txtDName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtDName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnDMake.
function btnDMake_Callback(hObject, eventdata, handles)
% hObject    handle to btnDMake (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    [ig1, realC] = implementC(handles);
    ctrlD = implementD(handles);
    
    h1 = figure;
    [ig1, ig2, w] = bode(ctrlD);     % Get proper frequency range
    bode(realC,w);
    grid;
    hold on;
    bode(ctrlD,w);

% --- Executes on button press in btnDSave.
function btnDSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnDSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    wsName = get(handles.txtDName, 'String');
    ctrlD = implementD(handles);
    assignin('base', wsName, ctrlD);



function txtCName_Callback(hObject, eventdata, handles)
% hObject    handle to txtCName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtCName as text
%        str2double(get(hObject,'String')) returns contents of txtCName as a double


% --- Executes during object creation, after setting all properties.
function txtCName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtCName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnCMake.
function btnCMake_Callback(hObject, eventdata, handles)
% hObject    handle to btnCMake (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    [ctrlC, realC] = implementC(handles);
    
    h1 = figure;
    [ig1, ig2, w] = bode(ctrlC);     % Get proper frequency range
    bode(realC, w);
    grid;
    hold on;
    bode(ctrlC,w);


% --- Executes on button press in btnCSave.
function btnCSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnCSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    wsName = get(handles.txtCName, 'String');
    ctrlC = implementC(handles);
    assignin('base', wsName, ctrlC);
    

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



function txtLam_Callback(hObject, eventdata, handles)
% hObject    handle to txtLam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtLam as text
%        str2double(get(hObject,'String')) returns contents of txtLam as a double


% --- Executes during object creation, after setting all properties.
function txtLam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtLam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtKd_Callback(hObject, eventdata, handles)
% hObject    handle to txtKd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtKd as text
%        str2double(get(hObject,'String')) returns contents of txtKd as a double


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



function txtWb_Callback(hObject, eventdata, handles)
% hObject    handle to txtWb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtWb as text
%        str2double(get(hObject,'String')) returns contents of txtWb as a double


% --- Executes during object creation, after setting all properties.
function txtWb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtWb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtWh_Callback(hObject, eventdata, handles)
% hObject    handle to txtWh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtWh as text
%        str2double(get(hObject,'String')) returns contents of txtWh as a double


% --- Executes during object creation, after setting all properties.
function txtWh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtWh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtN_Callback(hObject, eventdata, handles)
% hObject    handle to txtN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtN as text
%        str2double(get(hObject,'String')) returns contents of txtN as a double


% --- Executes during object creation, after setting all properties.
function txtN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menuType.
function menuType_Callback(hObject, eventdata, handles)
% hObject    handle to menuType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menuType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuType


% --- Executes during object creation, after setting all properties.
function menuType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menuCAppType.
function menuCAppType_Callback(hObject, eventdata, handles)
% hObject    handle to menuCAppType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menuCAppType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuCAppType


% --- Executes during object creation, after setting all properties.
function menuCAppType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuCAppType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menuDAppType.
function menuDAppType_Callback(hObject, eventdata, handles)
% hObject    handle to menuDAppType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menuDAppType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuDAppType


% --- Executes during object creation, after setting all properties.
function menuDAppType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuDAppType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menuC2D.
function menuC2D_Callback(hObject, eventdata, handles)
% hObject    handle to menuC2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menuC2D contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuC2D
if get(handles.menuC2D, 'Value') == 4
    set(handles.txtPrewarp, 'Enable', 'On');
else
    set(handles.txtPrewarp, 'Enable', 'Off');
end


% --- Executes during object creation, after setting all properties.
function menuC2D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuC2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtPrewarp_Callback(hObject, eventdata, handles)
% hObject    handle to txtPrewarp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtPrewarp as text
%        str2double(get(hObject,'String')) returns contents of txtPrewarp as a double


% --- Executes during object creation, after setting all properties.
function txtPrewarp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtPrewarp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtTs_Callback(hObject, eventdata, handles)
% hObject    handle to txtTs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTs as text
%        str2double(get(hObject,'String')) returns contents of txtTs as a double


% --- Executes during object creation, after setting all properties.
function txtTs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ----------------------------------------
% Main controller implementation functions
% ----------------------------------------
function [ctrlC, realC] = implementC(handles)

    % Get parameters
    Kp = str2num(get(handles.txtKp, 'String'));
    Ki = str2num(get(handles.txtKi, 'String'));
    Kd = str2num(get(handles.txtKd, 'String'));
    lam = str2num(get(handles.txtLam, 'String'));
    mu = str2num(get(handles.txtMu, 'String'));
    
    % Get real controller
    realC = fracpid(Kp,Ki,lam,Kd,mu);
    
    % Get approximation parameters
    wb = str2num(get(handles.txtWb, 'String'));
    wh = str2num(get(handles.txtWh, 'String'));
    N = str2num(get(handles.txtN, 'String'));
    allTypes = {'oust', 'ref'};
    type = allTypes{get(handles.menuType, 'Value')};
    
    % Get OustaPID
    ctrlC = oustapid(Kp,Ki,lam,Kd,mu,wb,wh,N,type);
    
    % Convert to desired model
    allModelTypes = {'zpk', 'tf', 'ss'};
    model = allModelTypes{get(handles.menuCAppType, 'Value')};
    
    % Return controller
    switch(model)
        case 'zpk'
            ctrlC = zpk(ctrlC);
        case 'tf'
            ctrlC = tf(ctrlC);
        case 'ss'
            ctrlC = ss(ctrlC);
    end

function ctrlD = implementD(handles)

    % Get continuous controller
    C1 = implementC(handles);
    
    % Get discretization method
    allC2D = {'zoh', 'foh', 'impulse', 'tustin', 'matched'};
    thisC2D = allC2D{get(handles.menuC2D, 'Value')};
    
    % Get prewarp freq
    prewarp = str2num(get(handles.txtPrewarp, 'String'));
    
    % Get sample time
    Ts = str2num(get(handles.txtTs, 'String'));
    
    % Obtain proper controller
    C1 = toproper(C1, str2num(get(handles.txtWh, 'String')));
    
    % Convert to discrete model
    if strcmp(thisC2D, 'tustin')
        % Use prewarp frequency
        opt = c2dOptions('Method', 'tustin', 'PrewarpFrequency', prewarp);
        ctrlD = c2d(C1, Ts, opt);
    else
        ctrlD = c2d(C1, Ts, thisC2D); 
    end
    
    % Convert to desired model
    allModelTypes = {'zpk', 'tf', 'ss'};
    model = allModelTypes{get(handles.menuDAppType, 'Value')};
    
    % Return controller
    switch(model)
        case 'zpk'
            ctrlD = zpk(ctrlD);
        case 'tf'
            ctrlD = tf(ctrlD);
        case 'ss'
            ctrlD = ss(ctrlD);
    end
    
    
