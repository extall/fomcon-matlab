function varargout = fid_pui(varargin)
% FID_PUI MATLAB code for fid_pui.fig
%      FID_PUI, by itself, creates a new FID_PUI or raises the existing
%      singleton*.
%
%      H = FID_PUI returns the handle to a new FID_PUI or the handle to
%      the existing singleton*.
%
%      FID_PUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FID_PUI.M with the given input arguments.
%
%      FID_PUI('Property','Value',...) creates a new FID_PUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fid_pui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fid_pui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fid_pui

% Last Modified by GUIDE v2.5 28-Jan-2016 18:37:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fid_pui_OpeningFcn, ...
                   'gui_OutputFcn',  @fid_pui_OutputFcn, ...
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


% --- Executes just before fid_pui is made visible.
function fid_pui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fid_pui (see VARARGIN)

% Choose default command line output for fid_pui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fid_pui wait for user response (see UIRESUME)
% uiwait(handles.figFidPui);

% Store the state of the buttons right away
% -1 means stopped, 0 means resume operation,
% 1 means STOP has been pressed, and 2 means Abort has been pressed
setappdata(handles.figFidPui, 'StopSignal', 0);
setappdata(handles.figFidPui, 'savedPerfData', []);


% --- Outputs from this function are returned to the command line.
function varargout = fid_pui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnStop.
function btnStop_Callback(hObject, eventdata, handles)
% hObject    handle to btnStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable both buttons and set the state
disableActions(hObject, eventdata, handles);
setStatusLine(hObject, eventdata, handles, ...
    'Identification stopped by user.');
setappdata(handles.figFidPui, 'StopSignal', 1);

% --- Executes on button press in btnTakeValues.
function btnTakeValues_Callback(hObject, eventdata, handles)
% hObject    handle to btnTakeValues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnResume.
function btnResume_Callback(hObject, eventdata, handles)
% hObject    handle to btnResume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function txtIterationNum_Callback(hObject, eventdata, handles)
% hObject    handle to txtIterationNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtIterationNum as text
%        str2double(get(hObject,'String')) returns contents of txtIterationNum as a double


% --- Executes during object creation, after setting all properties.
function txtIterationNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtIterationNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtCost_Callback(hObject, eventdata, handles)
% hObject    handle to txtCost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtCost as text
%        str2double(get(hObject,'String')) returns contents of txtCost as a double


% --- Executes during object creation, after setting all properties.
function txtCost_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtCost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtSimulationNum_Callback(hObject, eventdata, handles)
% hObject    handle to txtSimulationNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtSimulationNum as text
%        str2double(get(hObject,'String')) returns contents of txtSimulationNum as a double


% --- Executes during object creation, after setting all properties.
function txtSimulationNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSimulationNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtStepSize_Callback(hObject, eventdata, handles)
% hObject    handle to txtStepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtStepSize as text
%        str2double(get(hObject,'String')) returns contents of txtStepSize as a double


% --- Executes during object creation, after setting all properties.
function txtStepSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtStepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtFoOptimality_Callback(hObject, eventdata, handles)
% hObject    handle to txtFoOptimality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFoOptimality as text
%        str2double(get(hObject,'String')) returns contents of txtFoOptimality as a double


% --- Executes during object creation, after setting all properties.
function txtFoOptimality_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFoOptimality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtAlgorithm_Callback(hObject, eventdata, handles)
% hObject    handle to txtAlgorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAlgorithm as text
%        str2double(get(hObject,'String')) returns contents of txtAlgorithm as a double


% --- Executes during object creation, after setting all properties.
function txtAlgorithm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAlgorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menuView_Callback(hObject, eventdata, handles)
% hObject    handle to menuView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuPerformance_Callback(hObject, eventdata, handles)
% hObject    handle to menuPerformance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mySavedData = getappdata(handles.figFidPui, 'savedPerfData');
if ~isempty(mySavedData)
    h = figure;
    area(mySavedData.iterations, mySavedData.performance);
    xlabel('Iterations');
    ylabel('Performance (residual norm)');
    title('Identification performance graph');
else
    warndlg('There is no data to show!');
end


% --- Executes on button press in btnAbort.
function btnAbort_Callback(hObject, eventdata, handles)
% hObject    handle to btnAbort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable both buttons and set the state
disableActions(hObject, eventdata, handles);
setStatusLine(hObject, eventdata, handles, ...
   'Identification aborted by user.');
setappdata(handles.figFidPui, 'StopSignal', 2);

% --- USER DEFINED: Disables both action buttons
function disableActions(hObject, eventdata, handles)
set(handles.btnStop, 'Enable', 'off');
set(handles.btnAbort, 'Enable', 'off');

% --- USER DEFINED: Set the status line
function setStatusLine(hObject, eventdata, handles, text)
set(handles.txtIdentStatus, 'String', text); 


% --------------------------------------------------------------------
function menuFile_Callback(hObject, eventdata, handles)
% hObject    handle to menuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuSavePerfGraph_Callback(hObject, eventdata, handles)
% hObject    handle to menuSavePerfGraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Prompt user for workspace variable to export to
    
mySavedData = getappdata(handles.figFidPui, 'savedPerfData');
if ~isempty(mySavedData)
    name='Save Performance Graph to Workspace';
    prompt={'Workspace variable name:'};
    defaultanswer = {''};
    numlines=1;
    options.WindowStyle='normal';
    toExport=inputdlg(prompt,name,numlines,defaultanswer,options);
    
    if ~isempty(toExport)
        
        varName = toExport{1};
        assignin('base', varName, mySavedData)
        
    end
else
    warndlg('There is no data to save!');
end
