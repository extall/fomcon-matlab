function varargout = iopid_tune(varargin)
% IOPID_TUNE Integer-order PID tuning.
%
% Using this tool it is possible to obtain process models approximations
% from fractional-order models and to get corresponding controller gains
% by using analytical tuning formulae (e.g. the Ziegler-Nichols method).
%
% Last Modified by GUIDE v2.5 14-Feb-2011 00:46:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @iopid_tune_OpeningFcn, ...
                   'gui_OutputFcn',  @iopid_tune_OutputFcn, ...
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


% --- Executes just before iopid_tune is made visible.
function iopid_tune_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to iopid_tune (see VARARGIN)

% Fill in system name, if specified
sysName  = get(hObject, 'UserData');
if ~isempty(sysName)
    set(handles.txtPlant, 'String', sysName);
end

% Choose default command line output for iopid_tune
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes iopid_tune wait for user response (see UIRESUME)
% uiwait(handles.iopidTuningTool);


% --- Outputs from this function are returned to the command line.
function varargout = iopid_tune_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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


% --- Executes on selection change in menuMethod.
function menuMethod_Callback(hObject, eventdata, handles)
% hObject    handle to menuMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menuMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuMethod
    
    % Recompute gains
    ComputeGains(handles);


% --- Executes during object creation, after setting all properties.
function menuMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkPlotTune.
function chkPlotTune_Callback(hObject, eventdata, handles)
% hObject    handle to chkPlotTune (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkPlotTune


% --- Executes on button press in btnCompute.
function btnCompute_Callback(hObject, eventdata, handles)
% hObject    handle to btnCompute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    % Compute gains
    ComputeGains(handles);
    

% --- Executes on selection change in menuModel.
function menuModel_Callback(hObject, eventdata, handles)
% hObject    handle to menuModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menuModel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuModel
    
    modelType = get(hObject,'Value');
    switch modelType
        case 1
            % FOPDT
            tuneMethods = {'Ziegler-Nichols', ...
                           'Astrom-Hagglund (AMIGO)', ...
                           'Chien-Hrones-Reswick 1 (set-point regulation)', ...
                           'Chien-Hrones-Reswick 2 (disturbance rejection)', ...
                           'Cohen-Coon' };
                       
        case 2
            % IPDT
            tuneMethods = {'IPDT-ISE'};
            
        case 3
            % FOIPDT
            tuneMethods = {'FOIPDT tuning'};
    end
    
    set(handles.menuMethod, 'String', tuneMethods);
    set(handles.menuMethod, 'Value', 1);
    
    % Reset id fields
    set(handles.txtK, 'String', '50');
    set(handles.txtL, 'String', '50');
    set(handles.txtT, 'String', '50');
    
    % Reset tuned parameters
    set(handles.txtKp, 'String', '0');
    set(handles.txtKi, 'String', '0');
    set(handles.txtKd, 'String', '0');
    
    % Turn T off if model is 'IPDT'
    if modelType == 2
        set(handles.txtT, 'Enable', 'off');
    else
        set(handles.txtT, 'Enable', 'on');
    end


% --- Executes during object creation, after setting all properties.
function menuModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkPlotIdent.
function chkPlotIdent_Callback(hObject, eventdata, handles)
% hObject    handle to chkPlotIdent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkPlotIdent


% --- Executes on button press in btnIdentify.
function btnIdentify_Callback(hObject, eventdata, handles)
% hObject    handle to btnIdentify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    sysName = get(handles.txtPlant, 'String');
    [varex, varcl] = varexists(sysName);
    
    if varex && varcl
        
        % Plant model
        G = evalin('base', sysName);
        
        % Approximation
        appType = {'oust', 'ref'};
        appType = appType{get(handles.menuApprox,'Value')};
        
        w = str2num(get(handles.txtFreqRange,'String'));        
        N = str2num(get(handles.txtN,'String'));
        
        % Get model type
        models = {'fopdt', 'ipdt', 'foipdt'};
        modelString = models{get(handles.menuModel, 'Value')};
        
        % Get initial parameters
        K = str2num(get(handles.txtK,'String'));
        L = str2num(get(handles.txtL,'String'));
        T = str2num(get(handles.txtT,'String'));
        
        % Set simulation parameters
        sim = fsparam(G, appType, w, N);
        
        % Initial parameters
        initparams = [K L T];
        
        % Identify
        [K, L, T, y, t] = fotf2io(sim, modelString, initparams);
        
        % Write values
        set(handles.txtK, 'String', K);
        set(handles.txtL, 'String', L);
        set(handles.txtT, 'String', T);
        
        if get(handles.chkPlotIdent, 'Value')
           
            
            y_id = step(gpm(K,L,T,modelString), t);
            
            % Plot both initial system and the identified one
            h = figure();
            plot(t, y, t, y_id, '--', 'Linewidth', 2);
            set(h, 'NumberTitle', 'off');
            set(h, 'Name', [upper(modelString) ' identification results']);
            legend('Original model', ['Identified ' upper(modelString) ' model'], 'Location', 'Best');
            grid;
            
        end
        
    else
        
        errordlg('System no longer in workspace or invalid!');
        
    end


function txtPlant_Callback(hObject, eventdata, handles)
% hObject    handle to txtPlant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtPlant as text
%        str2double(get(hObject,'String')) returns contents of txtPlant as a double


% --- Executes during object creation, after setting all properties.
function txtPlant_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtPlant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtFreqRange_Callback(hObject, eventdata, handles)
% hObject    handle to txtFreqRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFreqRange as text
%        str2double(get(hObject,'String')) returns contents of txtFreqRange as a double


% --- Executes during object creation, after setting all properties.
function txtFreqRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFreqRange (see GCBO)
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


% --- Executes on selection change in menuApprox.
function menuApprox_Callback(hObject, eventdata, handles)
% hObject    handle to menuApprox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menuApprox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuApprox


% --- Executes during object creation, after setting all properties.
function menuApprox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuApprox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtK_Callback(hObject, eventdata, handles)
% hObject    handle to txtK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtK as text
%        str2double(get(hObject,'String')) returns contents of txtK as a double


% --- Executes during object creation, after setting all properties.
function txtK_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtL_Callback(hObject, eventdata, handles)
% hObject    handle to txtL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtL as text
%        str2double(get(hObject,'String')) returns contents of txtL as a double


% --- Executes during object creation, after setting all properties.
function txtL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtT_Callback(hObject, eventdata, handles)
% hObject    handle to txtT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtT as text
%        str2double(get(hObject,'String')) returns contents of txtT as a double


% --- Executes during object creation, after setting all properties.
function txtT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnTake.
function btnTake_Callback(hObject, eventdata, handles)
% hObject    handle to btnTake (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
% Open the FPID design tool
    ft = guidata(fpid('UserData', get(handles.txtPlant,'String')));
    
    % Set parameters (Kp,Ki,Kd)
    set(ft.txtKp,'String',get(handles.txtKp,'String'));
    set(ft.txtKi,'String',get(handles.txtKi,'String'));
    set(ft.txtKd,'String',get(handles.txtKd,'String'));
    
    % Set exponents
    set(ft.txtLambda,'String','1');
    set(ft.txtMu,'String','1');

function ComputeGains(handles)

    % Get parameters
    K = str2num(get(handles.txtK,'String'));
    L = str2num(get(handles.txtL,'String'));
    T = str2num(get(handles.txtT,'String'));
    
    switch get(handles.menuModel, 'Value')
        
        % Methods for FOPDT
        case 1
            
            switch get(handles.menuMethod, 'Value')

                case 1
                    % Ziegler-Nichols
                    [Kp,Ki,Kd] = zntune(K,L,T,3);

                case 2
                    % Astrom-Hagglund "AMIGO"
                    [Kp,Ki,Kd] = ahtune(K,L,T);

                case 3
                    % Chien-Hrones-Reswick 1 (set-point regulation)
                    [Kp,Ki,Kd] = chr1tune(K,L,T);

                case 4
                    % Chien-Hrones-Reswick 2 (disturbance rejection)
                    [Kp,Ki,Kd] = chr2tune(K,L,T);

                case 5
                    % Cohen-Coon
                    [Kp,Ki,Kd] = cctune(K,L,T);

                otherwise
                    errordlg('Unknown tuning method specified!', 'Error');

            end
        
        % Method for IPDT
        case 2
            
            [Kp,Ki,Kd] = ipdttune(K,L);
        
        % Method for FOIPDT
        case 3
            
            [Kp,Ki,Kd] = foipdttune(K,L,T);
    
    end
    
    % Set values
    set(handles.txtKp,'String',Kp);
    set(handles.txtKi,'String',Ki);
    set(handles.txtKd,'String',Kd);
