function varargout = fotfrid(varargin)
% FOTFRID Fractional-order transfer function frequency-domain identification.
%
% The tool allows to identify fractional-order transfer function from
% frequency response data.

% Last Modified by GUIDE v2.5 28-Apr-2011 21:51:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fotfrid_OpeningFcn, ...
                   'gui_OutputFcn',  @fotfrid_OutputFcn, ...
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


% --- Executes just before fotfrid is made visible.
function fotfrid_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fotfrid (see VARARGIN)

% Choose default command line output for fotfrid
handles.output = hObject;

% Initialize menu
refreshFfidataList(handles);

% Display first model
setIdModelGraphic(1, handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fotfrid wait for user response (see UIRESUME)
% uiwait(handles.guiFracIdent);


% --- Outputs from this function are returned to the command line.
function varargout = fotfrid_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function txtA_Callback(hObject, eventdata, handles)
% hObject    handle to txtA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtA as text
%        str2double(get(hObject,'String')) returns contents of txtA as a double


% --- Executes during object creation, after setting all properties.
function txtA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtA (see GCBO)
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


function txtB_Callback(hObject, eventdata, handles)
% hObject    handle to txtB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtB as text
%        str2double(get(hObject,'String')) returns contents of txtB as a double


% --- Executes during object creation, after setting all properties.
function txtB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnIdentify.
function btnIdentify_Callback(hObject, eventdata, handles)
% hObject    handle to btnIdentify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    freqIdentify(handles);

% --- Executes on button press in btnExport.
function btnExport_Callback(hObject, eventdata, handles)
% hObject    handle to btnExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Prompt user for workspace variable to export to
    name='Export Identified System to Workspace';
    prompt={'Workspace variable name:'};
    defaultanswer = {''};
    numlines=1;
    options.WindowStyle='normal';
    toExport=inputdlg(prompt,name,numlines,defaultanswer,options);
    
    if ~isempty(toExport)
        
        varName =toExport{1};
        G = newfotf(get(handles.txtB,'String'),get(handles.txtA,'String'));
        assignin('base', varName, G);
        
    end

function txtQ_Callback(hObject, eventdata, handles)
% hObject    handle to txtQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtQ as text
%        str2double(get(hObject,'String')) returns contents of txtQ as a double

% --- Executes during object creation, after setting all properties.
function txtQ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtM_Callback(hObject, eventdata, handles)
% hObject    handle to txtM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtM as text
%        str2double(get(hObject,'String')) returns contents of txtM as a double


% --- Executes during object creation, after setting all properties.
function txtM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in menuFidata.
function menuFidata_Callback(hObject, eventdata, handles)
% hObject    handle to menuFidata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menuFidata contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuFidata


% --- Executes during object creation, after setting all properties.
function menuFidata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuFidata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnRefresh.
function btnRefresh_Callback(hObject, eventdata, handles)
% hObject    handle to btnRefresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
refreshFfidataList(handles);


% --- Executes on selection change in menuPolySelect.
function menuPolySelect_Callback(hObject, eventdata, handles)
% hObject    handle to menuPolySelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menuPolySelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuPolySelect


% --- Executes during object creation, after setting all properties.
function menuPolySelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuPolySelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menuIdMethod.
function menuIdMethod_Callback(hObject, eventdata, handles)
% hObject    handle to menuIdMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menuIdMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuIdMethod


% --- Executes during object creation, after setting all properties.
function menuIdMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuIdMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Populate FIDATA list
function refreshFfidataList(handles)

    allVars = evalin('base', 'who');
    fdVars = {};
    
    for n=1:length(allVars)
        if strcmp(evalin('base',['class(' allVars{n} ')']),'ffidata')
            fdVars{end+1} = allVars{n};
        end
    end
    
    if ~isempty(fdVars)
       set(handles.menuFidata, 'String', fdVars);
       set(handles.menuFidata, 'Value', 1);
       set(handles.btnIdentify, 'Enable', 'on');
    else
       set(handles.menuFidata, 'String', ' ');
       set(handles.btnIdentify, 'Enable', 'off');
    end

function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to txtN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtN as text
%        str2double(get(hObject,'String')) returns contents of txtN as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menuIdType.
function menuIdType_Callback(hObject, eventdata, handles)
% hObject    handle to menuIdType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menuIdType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuIdType

    % Set corresponding model
    switch(get(hObject, 'Value'))
        case 1
            setIdModelGraphic(1, handles);
            set(handles.txtM, 'Enable', 'off');
            set(handles.menuTools, 'Enable', 'off');
        case {2, 3}
            setIdModelGraphic(2, handles);
            set(handles.txtM, 'Enable', 'on');
            set(handles.menuTools, 'Enable', 'on');
    end

% --- Executes during object creation, after setting all properties.
function menuIdType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuIdType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function setIdModelGraphic(modelNum, handles)

% Put image inside of picture axes
idModel = importdata(['ffid_mod' num2str(modelNum) '.jpg']);
axes(handles.axIdModel);
hi = imagesc(idModel);
axis off;

function freqIdentify(handles)

    % Get id data
    idd_all = get(handles.menuFidata,'String');
    idd_num = get(handles.menuFidata,'Value');
    idd = evalin('base', idd_all{idd_num});

    % Simulation method
    idmethods = {'h', 'l', 'v'};
    idmethod = idmethods{get(handles.menuIdType,'Value')};
    
    q = str2num(get(handles.txtQ, 'String'));
    
    n = str2num(get(handles.txtN, 'String'));
    m = str2num(get(handles.txtM, 'String'));
    
    % Identify model
    [a, na, b, nb, gg, J] = ffid(idd, q, [n, m], idmethod);
    
    % Update model
    G_zer = poly2str(b,nb);
    G_pol = poly2str(a,na);
    
    set(handles.txtB, 'String', G_zer);
    set(handles.txtA, 'String', G_pol);
    
    G_id = newfotf(G_zer, G_pol);
    
    [mag_id, ph_id] = bode(G_id, idd.w);
    
    % Squeeze
    mag_id = mag2db(squeeze(mag_id));
    ph_id = squeeze(ph_id);
    
	% Plot both initial system and the identified one
	h = figure();
    
    subplot(2,1,1);
    semilogx(idd.w, idd.mag, idd.w, mag_id, '--', 'Linewidth', 2);
    ylabel('Magnitude [dB]');
    legend('Initial data','Identified model','Location','Best');
    grid;
    
    subplot(2,1,2);
    semilogx(idd.w, idd.phase, idd.w, ph_id, '--', 'Linewidth', 2);
    xlabel('Frequency [rad/s]');
    ylabel('Phase [deg]');
    legend('Initial data','Identified model','Location','Best');
    grid;

    set(h, 'NumberTitle', 'off');
	set(h, 'Name', 'Frequency-domain identification results');
    
	config = fomcon('config');
    numSigDig = config.Core.General.Model_significant_digits;
	
    % Error index message
    msgbox(['Frequency-domain identification error index J = ' num2str(J, numSigDig)], 'Identification accuracy');

% --------------------------------------------------------------------
function menuTools_Callback(hObject, eventdata, handles)
% hObject    handle to menuTools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
           


% --------------------------------------------------------------------
function menuBestFit_Callback(hObject, eventdata, handles)
% hObject    handle to menuBestFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Input dialog
pr = inputdlg({'Max polynomial order:', 'Max iterations (0 - unlimited):'}, ...
               'Best fit optimization options', ...
               1, ...
               {'5', '0'});
               
if ~isempty(pr)
   
    % Get initial parameters
    idmethods = {'h', 'l', 'v'};
    idmethod = idmethods{get(handles.menuIdType,'Value')};
    
    % Get id data
    idd_all = get(handles.menuFidata,'String');
    idd_num = get(handles.menuFidata,'Value');
    idd = evalin('base', idd_all{idd_num});
    
    maxN = str2num(pr{1});
    maxiter = str2num(pr{2});
    
    % Initial guess
    qi = str2num(get(handles.txtQ, 'String'));
    ni = str2num(get(handles.txtN, 'String'));
    mi = str2num(get(handles.txtM, 'String'));
    
    init = [qi ni mi];
    
    % Pause to let dialog GUI terminate normally
    pause(2);
    
    % Optimize
    if maxiter ~= 0
        [q,n,m] = ffid_bf(idd, idmethod, init, maxN, maxiter);
    else
        [q,n,m] = ffid_bf(idd, idmethod, init, maxN);
    end
    
	config = fomcon('config');
    numSigDig = config.Core.General.Model_significant_digits;
	
    % Set new parameters
    set(handles.txtQ, 'String', num2str(q,numSigDig));
    set(handles.txtN, 'String', num2str(n,numSigDig));
    set(handles.txtM, 'String', num2str(m,numSigDig));
    
    % Run identification with these parameters
    freqIdentify(handles)
    
end
