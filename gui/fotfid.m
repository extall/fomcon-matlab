function varargout = fotfid(varargin)
% FOTFID Fractional-order transfer function identification tool.
%
% The tool is intended for identification of SISO systems by 
% fractional-order transfer function models.

% Last Modified by GUIDE v2.5 03-Feb-2016 23:11:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fotfid_OpeningFcn, ...
                   'gui_OutputFcn',  @fotfid_OutputFcn, ...
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


% --- Executes just before fotfid is made visible.
function fotfid_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fotfid (see VARARGIN)

% Choose default command line output for fotfid
handles.output = hObject;

% Initialize menu
refreshFidataList(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fotfid wait for user response (see UIRESUME)
% uiwait(handles.guiFracIdent);


% --- Outputs from this function are returned to the command line.
function varargout = fotfid_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on selection change in menuSimMethod.
function menuSimMethod_Callback(hObject, eventdata, handles)
% hObject    handle to menuSimMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menuSimMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuSimMethod
    
    % Toggle extra fields
    switch get(hObject, 'Value')
        case 1
            enEf = 'off';
        case {2, 3}
            enEf = 'on';
    end
    
    set(handles.txtFreqRange,'Enable',enEf);
    set(handles.txtN,'Enable',enEf);
    

% --- Executes during object creation, after setting all properties.
function menuSimMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuSimMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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
    
    % Get id data
    idd_all = get(handles.menuFidata,'String');
    idd_num = get(handles.menuFidata,'Value');
    idd = evalin('base', idd_all{idd_num});
    
    % Identification algorithm
    op = optimset();
    switch get(handles.menuAlgorithm, 'Value')
        case 1
            op.Display = 'iter';
        case 2
            op.Display = 'iter';
            op.IdentificationAlgorithm = 'lm';
            op.Lambda = str2num(get(handles.txtLambda,'String'));
    end
    
    % Initial system
    G = newfotf(get(handles.txtB,'String'),get(handles.txtA,'String'));
    
    % Delay parameter options
    if get(handles.chkL, 'Value')
       optL = str2num(get(handles.txtL, 'String'));
    else
       optL = [];
    end
    
    if (get(handles.chkK, 'Value'))
       optK = str2num(get(handles.txtK, 'String'));
    else
       optK = [];
    end

    % Simulation method
    simMethods = {'gl', 'oust', 'ref'};
    simMethod = simMethods{get(handles.menuSimMethod,'Value')};
    
    % Get parameters for approximations
    freq = str2num(get(handles.txtFreqRange,'String'));
    N = str2num(get(handles.txtN,'String'));
    
    switch(simMethod)
        case 'gl'
            fs = G;
        case 'oust'
            fs = fsparam(G, 'oust', freq, N);
        case 'ref'
            fs = fsparam(G, 'ref', freq, N);
    end
    
    % Coefficient limits
    lcb = str2num(get(handles.txtCMin,'String'));
    ucb = str2num(get(handles.txtCMax,'String'));
    
    % Exponent limits
    leb = str2num(get(handles.txtEMin,'String'));
    ueb = str2num(get(handles.txtEMax,'String'));
    
    % Limit coefficients?
    if get(handles.chkLTC,'Value')
        clim = [lcb; ucb];
    else
        clim = [];
    end
    
    % Limit exponents?
    if get(handles.chkLEX,'Value')
        elim = [leb; ueb];
    else
        elim = [];
    end
    
    % Polynomial fix
    polyfix = [get(handles.chkFixB,'Value');get(handles.chkFixA,'Value')];
    
    % Fix coefficients or exponents?
    identTypes = {'n', 'e', 'c'};
    identType = identTypes{get(handles.menuIdMethod, 'Value')};
    
    nPoints = [];
    pr = 0;
    
    % Ask about the number of points, if quick estimation is requested
    if get(handles.chkQuickEstimate, 'Value')
       pr = inputdlg({'Number of points to use for estimation:'}, ...
                   'Quick estimation', ...
                   1, ...
                   {num2str(floor(numel(idd.t)/2))});
               
       if ~isempty(pr)
           nPoints = str2num(pr{1});
       end
       
       % Pause to let dialog GUI terminate normally
       pause(2);
                   
    end
    
    if ~isempty(pr)
        [a, na, b, nb, L] = fid_(fs, {optK, optL}, idd, nPoints, ...
                                identType, polyfix, {clim, elim}, op);

        % Normalize displayed model
        if ~isempty(optK)
            K = dcgain(fotf(a,na,b,nb));
            b = b ./ K;
        end

        % Update model
        G_zer = poly2str(b,nb);
        G_pol = poly2str(a,na);

        set(handles.txtB, 'String', G_zer);
        set(handles.txtA, 'String', G_pol);

        G_id = newfotf(G_zer, G_pol);

		config = fomcon('config');
        numSigDig = config.Core.General.Model_significant_digits;
		
        if get(handles.chkK, 'Value')
            G_K = num2str(K, numSigDig);
            set(handles.txtK, 'String', G_K);
            G_id = G_id * str2num(G_K);
        end

        if get(handles.chkL, 'Value')
            G_L = num2str(L, numSigDig);
            G_id.ioDelay = L;
            set(handles.txtL, 'String', G_L);
        end

%
%      % Old behavior: use VALIDATE instead
%      %
%      %  y_id = lsim(G_id,idd.u,idd.t);
%      %
%      %  % Plot both initial system and the identified one
%      %  h = figure();
%      %  plot(idd.t, idd.y, idd.t, y_id, '--', 'Linewidth', 2);
%      %  hold on;
%      %  plot(idd.t, idd.y-y_id, 'Color', 'red');
%      %  hold off;
%      %  set(h, 'NumberTitle', 'off');
%      %  set(h, 'Name', 'Identification results');
%      %  legend('Source data', 'Identified model', 'Error', 'Location', 'Best');
%      %  grid;
    end

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
        
        varName = toExport{1};
        G = newfotf(get(handles.txtB,'String'),get(handles.txtA,'String'));
        % Add delay
        if (get(handles.chkL, 'Value'))
            G.ioDelay = str2num(get(handles.txtL, 'String'));
        end
        
        % Factor in static gain
        if (get(handles.chkK, 'Value'))
            G = G * str2num(get(handles.txtK, 'String'));
        end
        assignin('base', varName, G);
        
    end


function txtY_Callback(hObject, eventdata, handles)
% hObject    handle to txtY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtY as text
%        str2double(get(hObject,'String')) returns contents of txtY as a double


% --- Executes during object creation, after setting all properties.
function txtY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtY (see GCBO)
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



function txtModelN_Callback(hObject, eventdata, handles)
% hObject    handle to txtModelN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtModelN as text
%        str2double(get(hObject,'String')) returns contents of txtModelN as a double


% --- Executes during object creation, after setting all properties.
function txtModelN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtModelN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnGenerate.
function btnGenerate_Callback(hObject, eventdata, handles)
% hObject    handle to btnGenerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Check comm. order
    q = str2num(get(handles.txtQ,'String'));
    
    if q < 0.01 || q >= 2
        warndlg('Commensurate order out of acceptable range! Setting q=0.5 as default.', 'Error');
        q = 0.5;
        set(handles.txtQ,'String','0.5');
    end
    
    % Generate initial model
    N = str2num(get(handles.txtModelN, 'String'));
    
    na = q*(N:-1:1);
    na(end+1)=0;
    a = ones(1,length(na));
    
    % Get fractional polynomial
    poly = poly2str(a,na);
    
    % Set polynomial based on selection
    switch (get(handles.menuPolySelect, 'Value'))
        case 1
            set(handles.txtB,'String',poly);
        case 2
            set(handles.txtA,'String',poly);
        otherwise
            % Do nothing
    end
    
    drawnow();
    

% --- Executes on button press in chkFixExponents.
function chkFixExponents_Callback(hObject, eventdata, handles)
% hObject    handle to chkFixExponents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkFixExponents



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to txtB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtB as text
%        str2double(get(hObject,'String')) returns contents of txtB as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to txtA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtA as text
%        str2double(get(hObject,'String')) returns contents of txtA as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkLTC.
function chkLTC_Callback(hObject, eventdata, handles)
% hObject    handle to chkLTC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkLTC
if get(hObject,'Value')
   limEn = 'on'; 
else
   limEn = 'off';
end

set(handles.txtCMin,'Enable',limEn);
set(handles.txtCMax,'Enable',limEn);


function txtCMax_Callback(hObject, eventdata, handles)
% hObject    handle to txtCMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtCMax as text
%        str2double(get(hObject,'String')) returns contents of txtCMax as a double


% --- Executes during object creation, after setting all properties.
function txtCMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtCMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtCMin_Callback(hObject, eventdata, handles)
% hObject    handle to txtCMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtCMin as text
%        str2double(get(hObject,'String')) returns contents of txtCMin as a double


% --- Executes during object creation, after setting all properties.
function txtCMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtCMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtCAcc_Callback(hObject, eventdata, handles)
% hObject    handle to txtCAcc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtCAcc as text
%        str2double(get(hObject,'String')) returns contents of txtCAcc as a double


% --- Executes during object creation, after setting all properties.
function txtCAcc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtCAcc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtEAcc_Callback(hObject, eventdata, handles)
% hObject    handle to txtEAcc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtEAcc as text
%        str2double(get(hObject,'String')) returns contents of txtEAcc as a double


% --- Executes during object creation, after setting all properties.
function txtEAcc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtEAcc (see GCBO)
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


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
refreshFidataList(handles);


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


% --- Executes on button press in chkLEX.
function chkLEX_Callback(hObject, eventdata, handles)
% hObject    handle to chkLEX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkLEX
switch(get(hObject, 'Value'))
    case 0
        set(handles.txtEMin,'Enable','off');
        set(handles.txtEMax,'Enable','off');

    case 1
        set(handles.txtEMin,'Enable','on');
        set(handles.txtEMax,'Enable','on');
end


function txtEMax_Callback(hObject, eventdata, handles)
% hObject    handle to txtEMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtEMax as text
%        str2double(get(hObject,'String')) returns contents of txtEMax as a double


% --- Executes during object creation, after setting all properties.
function txtEMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtEMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtEMin_Callback(hObject, eventdata, handles)
% hObject    handle to txtEMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtEMin as text
%        str2double(get(hObject,'String')) returns contents of txtEMin as a double


% --- Executes during object creation, after setting all properties.
function txtEMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtEMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkFixA.
function chkFixA_Callback(hObject, eventdata, handles)
% hObject    handle to chkFixA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkFixA
switch get(hObject,'Value')
    case 0
        set(handles.txtA, 'Enable', 'on');
    case 1
        set(handles.txtA, 'Enable', 'off');
        
        % Check the other control
        if get(handles.chkFixB, 'Value')
            set(handles.chkFixB, 'Value', 0);
            set(handles.txtB,'Enable','on');
        end
end

checkGainSetting(handles)

% --- Executes on button press in chkFixB.
function chkFixB_Callback(hObject, eventdata, handles)
% hObject    handle to chkFixB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkFixB
switch get(hObject,'Value')
    case 0
        set(handles.txtB, 'Enable', 'on');
    case 1
        set(handles.txtB, 'Enable', 'off');
        
        % Check the other control
        if get(handles.chkFixA, 'Value')
            set(handles.chkFixA, 'Value', 0);
            set(handles.txtA,'Enable','on');
        end
end

checkGainSetting(handles)


% --- Executes on selection change in menuIdMethod.
function menuIdMethod_Callback(hObject, eventdata, handles)
% hObject    handle to menuIdMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menuIdMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuIdMethod

checkGainSetting(handles)


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
function refreshFidataList(handles)

    allVars = evalin('base', 'who');
    fdVars = {};
    
    for n=1:length(allVars)
        if strcmp(evalin('base',['class(' allVars{n} ')']),'fidata')
            fdVars{end+1} = allVars{n};
        end
    end
    
    if ~isempty(fdVars)
       set(handles.menuFidata, 'String', fdVars);
       set(handles.menuFidata, 'Value', 1);
       set(handles.btnIdentify, 'Enable', 'on');
       set(handles.menuData, 'Enable', 'on');
    else
       set(handles.menuFidata, 'String', ' ');
       set(handles.btnIdentify, 'Enable', 'off');
       set(handles.menuData, 'Enable', 'off');
    end
    
    
% Check gain manual setting
function checkGainSetting(handles)
    
    if ~(get(handles.chkFixB, 'Value') == 1) && ...
       ~(get(handles.chkFixA, 'Value') == 1) && ...
        (get(handles.menuIdMethod, 'Value') == 1)
            set(handles.chkK, 'Enable', 'on');
    else
        set(handles.chkK, 'Enable', 'off');
        set(handles.chkK, 'Value', 0);
        set(handles.txtK, 'Enable', 'off');
    end

% --------------------------------------------------------------------
function menuImport_Callback(hObject, eventdata, handles)
% hObject    handle to menuImport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuIGModel_Callback(hObject, eventdata, handles)
% hObject    handle to menuIGModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pr = inputdlg({'Model workspace name:'}, ...
               'Import initial guess model', ...
               1, ...
               {''});
           
if ~isempty(pr)
   modelName = pr{1};
   myModel = evalin('base', modelName);
   [a,na,b,nb] = fotfparam(myModel);
   set(handles.txtB, 'String', poly2str(b,nb));
   set(handles.txtA, 'String', poly2str(a,na));
   
   % Disable DC gain setting : TODO: possibility to check model for DC gain
   set(handles.chkK, 'Value', 0);
   chkK_Callback(handles.chkK, [], handles);
end


% --- Executes on button press in chkQuickEstimate.
function chkQuickEstimate_Callback(hObject, eventdata, handles)
% hObject    handle to chkQuickEstimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkQuickEstimate


% --- Executes on button press in chkK.
function chkK_Callback(hObject, eventdata, handles)
% hObject    handle to chkK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkK

% Enable/disable textbox
switch(get(handles.chkK, 'Value'))
    case 0
        set(handles.txtK, 'Enable', 'Off');
    case 1
        set(handles.txtK, 'Enable', 'On');
    otherwise
        % What can be done otherwise?
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


% --- Executes on button press in chkL.
function chkL_Callback(hObject, eventdata, handles)
% hObject    handle to chkL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkL

% Enable/disable textbox
switch(get(handles.chkL, 'Value'))
    case 0
        set(handles.txtL, 'Enable', 'Off');
    case 1
        set(handles.txtL, 'Enable', 'On');
    otherwise
        % What can be done otherwise?
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

% --------------------------------------------------------------------
function menuImportTD_Callback(hObject, eventdata, handles)
% hObject    handle to menuImportTD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pr = inputdlg({'Output vector:', 'Input vector:', 'Time vector:', 'New data object name:'}, ...
               'Import time domain data from workspace', ...
               1, ...
               {'y','u','t','id_'});
           
if ~isempty(pr)
   y = evalin('base', pr{1});
   u = evalin('base', pr{2});
   t = evalin('base', pr{3});
   assignin('base', pr{4}, fidata(y,u,t));
   refreshFidataList(handles);
   
   allDatasets = get(handles.menuFidata, 'String');
   for i=1:length(allDatasets)
      if strcmpi(allDatasets{i},pr{4})
          set(handles.menuFidata, 'Value', i);
      end
   end
end


% --------------------------------------------------------------------
function menuData_Callback(hObject, eventdata, handles)
% hObject    handle to menuData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuPlot_Callback(hObject, eventdata, handles)
% hObject    handle to menuPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get id data
idd_all = get(handles.menuFidata,'String');
idd_num = get(handles.menuFidata,'Value');
idd = evalin('base', idd_all{idd_num});
    
% Plot the data
plot(idd);

% --------------------------------------------------------------------
function menuValidate_Callback(hObject, eventdata, handles)
% hObject    handle to menuValidate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get id data
idd_all = get(handles.menuFidata,'String');
idd_num = get(handles.menuFidata,'Value');
idd = evalin('base', idd_all{idd_num});

% Generate the model
G = newfotf(get(handles.txtB,'String'),get(handles.txtA,'String'));

% Add delay
if (get(handles.chkL, 'Value'))
   G.ioDelay = str2num(get(handles.txtL, 'String'));
end
        
% Factor in static gain
if (get(handles.chkK, 'Value'))
   G = G * str2num(get(handles.txtK, 'String'));
end

% Simulation method
simMethods = {'gl', 'oust', 'ref'};
simMethod = simMethods{get(handles.menuSimMethod,'Value')};

% Get parameters for approximations
freq = str2num(get(handles.txtFreqRange,'String'));
N = str2num(get(handles.txtN,'String'));

switch(simMethod)
    case 'gl'
        fs = G;
    case 'oust'
        fs = fsparam(G, 'oust', freq, N);
    case 'ref'
        fs = fsparam(G, 'ref', freq, N);
end

% Plot the validation graph
validate(idd,fs);


% --------------------------------------------------------------------
function menuModel_Callback(hObject, eventdata, handles)
% hObject    handle to menuModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuNormalize_Callback(hObject, eventdata, handles)
% hObject    handle to menuNormalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.chkK, 'Value')
    K = str2num(get(handles.txtK, 'String'));
else
    K = 1;
end

% Normalize model
G = normalize(K*newfotf(get(handles.txtB,'String'),get(handles.txtA,'String')));
[a,na,b,nb]=fotfparam(G);

if get(handles.chkK, 'Value')
   dcGain = dcgain(G);
   if ~fleq(dcGain,0) && ~isnan(dcGain) && ~isinf(dcGain)
       set(handles.txtK, 'String', dcgain(G));
       b = b ./ dcgain(G);
   else
       set(handles.txtK, 'String', '1');
   end
end

set(handles.txtB,'String',poly2str(b,nb));
set(handles.txtA,'String',poly2str(a,na));


% --------------------------------------------------------------------
function menuRound_Callback(hObject, eventdata, handles)
% hObject    handle to menuRound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pr = inputdlg({'Coefficient accuracy:', 'Exponent accuracy:'}, ...
               'Round model parameters', ...
               1, ...
               {'0.001','0.001'});
           
if ~isempty(pr)
    G=newfotf(get(handles.txtB,'String'),get(handles.txtA,'String'));
    [a,na,b,nb]=fotfparam(round(G,str2num(pr{2}),str2num(pr{1})));
    set(handles.txtB,'String',poly2str(b,nb));
    set(handles.txtA,'String',poly2str(a,na));
end


% --------------------------------------------------------------------
function menuTruncate_Callback(hObject, eventdata, handles)
% hObject    handle to menuTruncate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pr = inputdlg({'Coefficient accuracy:', 'Exponent accuracy:'}, ...
               'Truncate model parameters', ...
               1, ...
               {'0.001','0.001'});
           
if ~isempty(pr)
    G=newfotf(get(handles.txtB,'String'),get(handles.txtA,'String'));
    [a,na,b,nb]=fotfparam(trunc(G,str2num(pr{2}),str2num(pr{1})));
    set(handles.txtB,'String',poly2str(b,nb));
    set(handles.txtA,'String',poly2str(a,na));
end


% --------------------------------------------------------------------
function menuIsStable_Callback(hObject, eventdata, handles)
% hObject    handle to menuIsStable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Generate the model
G = newfotf(get(handles.txtB,'String'),get(handles.txtA,'String'));

% Add delay
if (get(handles.chkL, 'Value'))
   G.ioDelay = str2num(get(handles.txtL, 'String'));
end
        
% Factor in static gain
if (get(handles.chkK, 'Value'))
   G = G * str2num(get(handles.txtK, 'String'));
end

config = fomcon('config');
numSigDig = config.Core.General.Model_significant_digits;

% Create figure
h  = figure;
ax = axes;

% Check stability
[K,q,err,apol] = isstable(G, true);

config = fomcon('config');
numSigDig = config.Core.General.Model_significant_digits;

% Display message box with resolution
if (K == 1)
    msgbox(['System appears to be STABLE with order q=' num2str(q,numSigDig)], ...
        ['Stability test for '' G '''], 'modal');
else
    msgbox(['System appears to be UNSTABLE with order q=' num2str(q,numSigDig)], ...
        ['Stability test for '' G '''], 'modal');
end

% Set figure name
set(h, 'NumberTitle', 'off');
set(h, 'Name', [''' G ''' ' ' 'Stability Test Results']);


% --------------------------------------------------------------------
function menuResiduals_Callback(hObject, eventdata, handles)
% hObject    handle to menuResiduals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get id data
idd_all = get(handles.menuFidata,'String');
idd_num = get(handles.menuFidata,'Value');
idd = evalin('base', idd_all{idd_num});

% Generate the model
G = newfotf(get(handles.txtB,'String'),get(handles.txtA,'String'));

% Add delay
if (get(handles.chkL, 'Value'))
   G.ioDelay = str2num(get(handles.txtL, 'String'));
end
        
% Factor in static gain
if (get(handles.chkK, 'Value'))
   G = G * str2num(get(handles.txtK, 'String'));
end

% Simulation method
simMethods = {'gl', 'oust', 'ref'};
simMethod = simMethods{get(handles.menuSimMethod,'Value')};

% Get parameters for approximations
freq = str2num(get(handles.txtFreqRange,'String'));
N = str2num(get(handles.txtN,'String'));

switch(simMethod)
    case 'gl'
        fs = G;
    case 'oust'
        fs = fsparam(G, 'oust', freq, N);
    case 'ref'
        fs = fsparam(G, 'ref', freq, N);
end

% Plot the validation graph
resids(idd,fs);


% --- Executes on selection change in menuAlgorithm.
function menuAlgorithm_Callback(hObject, eventdata, handles)
% hObject    handle to menuAlgorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menuAlgorithm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuAlgorithm
alg = get(handles.menuAlgorithm, 'Value');
if (alg==1)
    
    set(handles.txtLambda, 'Enable', 'off');
    
    % Enable bounds
    set(handles.chkLTC, 'Enable', 'on');
    set(handles.chkLEX, 'Enable', 'on');
    
    % If corresponding boxes checked, enable them
    if (get(handles.chkLTC, 'Value'))
        set(handles.txtCMin, 'Enable', 'on');
        set(handles.txtCMax, 'Enable', 'on');
    end
    
    if (get(handles.chkLEX, 'Value'))
        set(handles.txtEMin, 'Enable', 'on');
        set(handles.txtEMax, 'Enable', 'on');
    end
    
else
    
    set(handles.txtLambda, 'Enable', 'on');
    
    % Disable bounds
    set(handles.chkLTC, 'Enable', 'off');
    set(handles.chkLEX, 'Enable', 'off');
    
    % Disable all other boxes as well
    set(handles.txtCMin, 'Enable', 'off');
    set(handles.txtCMax, 'Enable', 'off');
    set(handles.txtEMin, 'Enable', 'off');
    set(handles.txtEMax, 'Enable', 'off');
    
end


% --- Executes during object creation, after setting all properties.
function menuAlgorithm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuAlgorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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


% --------------------------------------------------------------------
function menuUncertaintyTest_Callback(hObject, eventdata, handles)
% hObject    handle to menuUncertaintyTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pr = inputdlg({'Method (''oat'' or ''monte-carlo'')', ...
               'Coefficient uncertainty {0.0, 1.0}', ...
               'Exponent uncertainty {0.0, 1.0}', ...
               'Simulations to perform (for ''oat'' only)'}, ...
               'Uncertainty test', ...
               1, ...
               {'oat', '0.01', '0.01', '25'});
           
if ~isempty(pr)

    % Get id data
    idd_all = get(handles.menuFidata,'String');
    idd_num = get(handles.menuFidata,'Value');
    idd = evalin('base', idd_all{idd_num});
    
    % Generate the model
    G = newfotf(get(handles.txtB,'String'),get(handles.txtA,'String'));
    
    % Add delay
    if (get(handles.chkL, 'Value'))
        G.ioDelay = str2num(get(handles.txtL, 'String'));
    end
    
    % Factor in static gain
    if (get(handles.chkK, 'Value'))
        G = G * str2num(get(handles.txtK, 'String'));
    end
    
    % Simulation method
    simMethods = {'gl', 'oust', 'ref'};
    simMethod = simMethods{get(handles.menuSimMethod,'Value')};
    
    % Get parameters for approximations
    freq = str2num(get(handles.txtFreqRange,'String'));
    N = str2num(get(handles.txtN,'String'));
    
    switch(simMethod)
        case 'gl'
            fs = G;
        case 'oust'
            fs = fsparam(G, 'oust', freq, N);
        case 'ref'
            fs = fsparam(G, 'ref', freq, N);
    end
    
    % Call the testsim function
    testsim(idd,fs,pr{1},[str2double(pr{2}), str2double(pr{3})], ...
        str2double(pr{4}));
    
end


% --------------------------------------------------------------------
function menuModelActions_Callback(hObject, eventdata, handles)
% hObject    handle to menuModelActions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
