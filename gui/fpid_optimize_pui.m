function varargout = fpid_optimize_pui(varargin)
% FPID_OPTIMIZE_PUI MATLAB code for fpid_optimize_pui.fig
%      FPID_OPTIMIZE_PUI, by itself, creates a new FPID_OPTIMIZE_PUI or raises the existing
%      singleton*.
%
%      H = FPID_OPTIMIZE_PUI returns the handle to a new FPID_OPTIMIZE_PUI or the handle to
%      the existing singleton*.
%
%      FPID_OPTIMIZE_PUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FPID_OPTIMIZE_PUI.M with the given input arguments.
%
%      FPID_OPTIMIZE_PUI('Property','Value',...) creates a new FPID_OPTIMIZE_PUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fpid_optimize_pui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fpid_optimize_pui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fpid_optimize_pui

% Last Modified by GUIDE v2.5 02-Feb-2016 23:05:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fpid_optimize_pui_OpeningFcn, ...
                   'gui_OutputFcn',  @fpid_optimize_pui_OutputFcn, ...
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


% --- Executes just before fpid_optimize_pui is made visible.
function fpid_optimize_pui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fpid_optimize_pui (see VARARGIN)

% Choose default command line output for fpid_optimize_pui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fpid_optimize_pui wait for user response (see UIRESUME)
% uiwait(handles.figFpidOptimizePui);

% Store the state of the buttons right away
% -1 means stopped, 0 means resume operation,
% 1 means STOP has been pressed, and 2 means Abort has been pressed
setappdata(handles.figFpidOptimizePui, 'StopSignal', 0);
setappdata(handles.figFpidOptimizePui, 'savedPerfData', []);


% --- Outputs from this function are returned to the command line.
function varargout = fpid_optimize_pui_OutputFcn(hObject, eventdata, handles) 
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
    'Optimization stopped by user.');
setappdata(handles.figFpidOptimizePui, 'StopSignal', 1);

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



function txtPerformance_Callback(hObject, eventdata, handles)
% hObject    handle to txtPerformance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtPerformance as text
%        str2double(get(hObject,'String')) returns contents of txtPerformance as a double


% --- Executes during object creation, after setting all properties.
function txtPerformance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtPerformance (see GCBO)
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



function txtOptimalityType_Callback(hObject, eventdata, handles)
% hObject    handle to txtOptimalityType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtOptimalityType as text
%        str2double(get(hObject,'String')) returns contents of txtOptimalityType as a double


% --- Executes during object creation, after setting all properties.
function txtOptimalityType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtOptimalityType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtStepNorm_Callback(hObject, eventdata, handles)
% hObject    handle to txtStepNorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtStepNorm as text
%        str2double(get(hObject,'String')) returns contents of txtStepNorm as a double


% --- Executes during object creation, after setting all properties.
function txtStepNorm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtStepNorm (see GCBO)
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
function menuViewPerformance_Callback(hObject, eventdata, handles)
% hObject    handle to menuViewPerformance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mySavedData = getappdata(handles.figFpidOptimizePui, 'savedPerfData');
if ~isempty(mySavedData)
    h = figure;
    title('Optimization performance dynamics');
    
    % Draw the graphs for parameter dynamics
    subplot(3,1,1);
    plot(mySavedData.iterations, mySavedData.parameters(:,1),'-r*');
    hold on;
    plot(mySavedData.iterations, mySavedData.parameters(:,2),'-.bo');
    plot(mySavedData.iterations, mySavedData.parameters(:,3),'--ks');
    legend('Kp','Ki','Kd');
    ylabel('Gains');
    xlim([min(mySavedData.iterations) max(mySavedData.iterations)]);
    grid;
    
    subplot(3,1,2);
    plot(mySavedData.iterations, mySavedData.parameters(:,4),'-.bo');
    hold on;
    plot(mySavedData.iterations, mySavedData.parameters(:,5),'--ks');
    legend('\lambda','\mu');
    ylabel('Orders');
    xlim([min(mySavedData.iterations) max(mySavedData.iterations)]);
    grid;
    
    subplot(3,1,3);
    area(mySavedData.iterations, mySavedData.performance);
    ylabel(['Performance: ' mySavedData.metric]);
    xlabel('Iterations');
    xlim([min(mySavedData.iterations) max(mySavedData.iterations)]);
    
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
   'Optimization aborted by user.');
setappdata(handles.figFpidOptimizePui, 'StopSignal', 2);

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
    
mySavedData = getappdata(handles.figFpidOptimizePui, 'savedPerfData');
if ~isempty(mySavedData)
    name='Save Performance Dynamics to Workspace';
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



function txtFeasibility_Callback(hObject, eventdata, handles)
% hObject    handle to txtFeasibility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFeasibility as text
%        str2double(get(hObject,'String')) returns contents of txtFeasibility as a double


% --- Executes during object creation, after setting all properties.
function txtFeasibility_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFeasibility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
