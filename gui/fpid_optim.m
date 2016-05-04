function varargout = fpid_optim(varargin)
% FPID_OPTIM Fractional-order PID optimization tool.
%
% The tool allows to obtain a set of optimized controller parameters given
% a set of design specifications.
%
% Last Modified by GUIDE v2.5 05-Nov-2013 16:34:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fpid_optim_OpeningFcn, ...
                   'gui_OutputFcn',  @fpid_optim_OutputFcn, ...
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

% --- Executes just before fpid_optim is made visible.
function fpid_optim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fpid_optim (see VARARGIN)

% Get optimization data if supplied
toOpt  = get(hObject, 'UserData');

if ~isempty(toOpt)
    set(handles.txtPlant, 'String', toOpt.G);
    set(handles.txtKp, 'String', toOpt.Kp);
    set(handles.txtKi, 'String', toOpt.Ki);
    set(handles.txtKd, 'String', toOpt.Kd);
    set(handles.txtLam, 'String', toOpt.lam);
    set(handles.txtMu, 'String', toOpt.mu);
end

% Initialize UserData
if isempty(get(hObject, 'UserData'))
	set(hObject, 'UserData', []);
end

% Check plant type
checkPlantType(handles);

% Choose default command line output for fpid_optim
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fpid_optim wait for user response (see UIRESUME)
% uiwait(handles.figFpidOptimTool);


% --- Outputs from this function are returned to the command line.
function varargout = fpid_optim_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in menuMetric.
function menuMetric_Callback(hObject, eventdata, handles)
% hObject    handle to menuMetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menuMetric contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuMetric


% --- Executes during object creation, after setting all properties.
function menuMetric_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuMetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnOptimize.
function btnOptimize_Callback(hObject, eventdata, handles)
% hObject    handle to btnOptimize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    checkPlantType(handles);

    sysName = get(handles.txtPlant, 'String');
    [varex, varcl, varclass] = varexists(sysName);
    
    if varex
       
        % Get all parameters
        
        % Check system type
        switch(lower(varclass))
            case {'tf', 'zpk', 'ss'}
                nonfotf = true;
            case 'fotf'
                nonfotf = false;
            otherwise
                errordlg('Unknown LTI model type!', 'Error');
        end

        % Plant model
        G = evalin('base', sysName);
        
        % Approximation
        appType = {'oust', 'ref'};
        type = appType{get(handles.menuApprox,'Value')};
        
        % Convert objects to state-space
        if strcmpi(varclass, 'tf') || strcmpi(varclass, 'zpk')
            G = ss(G);
        end
        
        w = str2num(get(handles.txtFreqRange,'String'));
        
        N = str2num(get(handles.txtAppOrder,'String'));
        
		config = fomcon('config');
        numSigDig = config.Core.General.Model_significant_digits;
		
        % Get FPID parameters
        Kp_val = get(handles.txtKp,'String');
        if ~isempty(str2num(Kp_val))
            Kp.val = str2num(get(handles.txtKp,'String'));
        else
            Kp.val = evalin('base', Kp_val);
            set(handles.txtKp,'String',num2str(Kp.val,numSigDig));
        end
        Kp.min = str2num(get(handles.txtKpMin,'String'));
        Kp.max = str2num(get(handles.txtKpMax,'String'));
        
        Ki_val = get(handles.txtKi,'String');
        if ~isempty(str2num(Ki_val))
            Ki.val = str2num(get(handles.txtKi,'String'));
        else
            Ki.val = evalin('base', Ki_val);
            set(handles.txtKi,'String',num2str(Ki.val,numSigDig));
        end
        Ki.min = str2num(get(handles.txtKiMin,'String'));
        Ki.max = str2num(get(handles.txtKiMax,'String'));
        
        Kd_val = get(handles.txtKd,'String');
        if ~isempty(str2num(Kd_val))
            Kd.val = str2num(get(handles.txtKd,'String'));
        else
            Kd.val = evalin('base', Kd_val);
            set(handles.txtKd,'String',num2str(Kd.val,numSigDig));
        end
        Kd.min = str2num(get(handles.txtKdMin,'String'));
        Kd.max = str2num(get(handles.txtKdMax,'String'));
        
        % Fractional exponents
        lam_val = get(handles.txtLam,'String');
        if ~isempty(str2num(lam_val))
            lam.val = str2num(get(handles.txtLam,'String'));
        else
            lam.val = evalin('base', lam_val);
            set(handles.txtLam,'String',num2str(lam.val,numSigDig));
        end
        lam.min = str2num(get(handles.txtLamMin,'String'));
        lam.max = str2num(get(handles.txtLamMax,'String'));
        
        mu_val = get(handles.txtMu,'String');
        if ~isempty(str2num(mu_val))
            mu.val = str2num(get(handles.txtMu,'String'));
        else
            mu.val = evalin('base', mu_val);
            set(handles.txtMu,'String',num2str(mu.val,numSigDig));
        end
        mu.min = str2num(get(handles.txtMuMin,'String'));
        mu.max = str2num(get(handles.txtMuMax,'String'));
        
        % Get performance metric
        per    = {'ise','iae','itse','itae'};
        metric = per{get(handles.menuMetric,'Value')};
        
        % Get algorithm
        allAlgs = {'nelder-mead', 'interior-point', 'sqp', 'active-set'};
        alg     = allAlgs{get(handles.menuOptimizerAlgorithm, 'Value')};
        
        % Get constraints
        cancelzero = get(handles.chkZeroCancelation, 'Value');
        
        % Maximum simulation time
        maxtime = str2num(get(handles.txtMaxTime, 'String'));
        dtmin   = str2num(get(handles.txtDtMin, 'String'));
        dtmax   = str2num(get(handles.txtDtMax, 'String'));
        
        simtime = [dtmin; dtmax; maxtime];
            
        % Margins
        if get(handles.chkGMEnable, 'Value')
            margins = [str2num(get(handles.txtGm, 'String')), ...
                       get(handles.chkGmExact, 'Value');
                       str2num(get(handles.txtPm, 'String')), ...
                       get(handles.chkPmExact, 'Value')];
        else
            margins = [];
        end
        
        % Sensitivity
        if get(handles.chkSensitivity, 'Value')
            sens = [str2num(get(handles.txtWt,  'String')), ...
                    str2num(get(handles.txtTdb, 'String'));
                    str2num(get(handles.txtWs,  'String')), ...
                    str2num(get(handles.txtSdb, 'String'))];
        else
            sens = [];
        end
        
        % Gain variations
        if get(handles.chkEnableRobustness, 'Value')
            if ~get(handles.chkWh, 'Value')
                gainvar = [str2num(get(handles.txtW, 'String'));
                           str2num(get(handles.txtW, 'String'))]; 
            else
                gainvar = [str2num(get(handles.txtW, 'String'));
                           str2num(get(handles.txtRange, 'String'))];
            end
        else
            gainvar = [];
        end
        
        % Control signal limits
        if get(handles.chkControlLimits, 'Value')
            ulim = [str2num(get(handles.txtUMin, 'String'));
                    str2num(get(handles.txtUMax, 'String'))];
        else
            ulim = [];
        end
        
        % Control law weighting
        if get(handles.chkWeight, 'Value')
            wgt = str2num(get(handles.txtWeight, 'String'));
        else
            wgt = [];
        end
        
        % Setpoint
        sp = evalin('base',get(handles.txtSP, 'String'));
        
        % Check whether "use Simulink" is enabled
        % Enable, if SP is a structure
        if (isa(sp, 'struct') && ~get(handles.chkUseSimulink, 'Value'))
           errordlg('Using a custom reference signal also requires using Simulink for system simulation');
           return;
        end
        
        % Strict constraints?
        stct =  get(handles.chkStrict,'Value');
        
        % Optimize PID control
        op = optimset('display','iter');
        
        % Check if MaxIter is to be set
        if get(handles.chkLimitIter,'Value')
           op.MaxIter = str2num(get(handles.txtLimitIter,'String')); 
        end
        
        optimType = {'n', 'e', 'g'};
        optType = optimType{get(handles.menuFix, 'Value')};
        
        % Build parameter structures
        if nonfotf
            fsim = fsparam(fotf(), type, w, N);
        else
            fsim = fsparam(G, type, w, N);
        end

        % Optimization options
        fopt = fpopt(optType, ...
                     [Kp.val Ki.val Kd.val lam.val mu.val], ...
                     [Kp.max Ki.max Kd.max lam.max mu.max], ...
                     [Kp.min Ki.min Kd.min lam.min mu.min], ...
                     metric, alg, sp, margins, sens, ulim,  ...
                     wgt, gainvar, simtime, cancelzero, ...
                     stct, op);

        % Use Simulink for simulation?
        useSimChk = get(handles.chkUseSimulink, 'Value');
        
        if useSimChk
            allModels = get(handles.menuSimulinkModel, 'String');
            usesim = allModels{get(handles.menuSimulinkModel, 'Value')};
        else
            usesim = [];
        end
            
        % Pause to let dialog GUI terminate normally
       pause(2);
        
        % Obtain optimized parameters
        %try
            if nonfotf
                [Kp, Ki, Kd, lam, mu, Z] = fpid_optim_fpid_optimize(fsim, fopt, G, usesim, handles);
            else     
                [Kp, Ki, Kd, lam ,mu, Z] = fpid_optim_fpid_optimize(fsim, fopt, [], usesim, handles);
            end
        %catch e
        %    errordlg(['An error occured during optimization:' e.message],'Error');
        %    return;
        %end
        
        % Set new parameters
        set(handles.txtKp, 'String', num2str(Kp,numSigDig));
        set(handles.txtKi, 'String', num2str(Ki,numSigDig));
        set(handles.txtKd, 'String', num2str(Kd,numSigDig));
        set(handles.txtLam, 'String', num2str(lam,numSigDig));
        set(handles.txtMu, 'String', num2str(mu,numSigDig));
        
        if get(handles.chkSimulate,'Value') || ...
           get(handles.chkSimulateOnly,'Value')
           
           t = Z.timeSpan;
           
           if ~nonfotf
               G = oustapp(G, w(1), w(2), N, type);
           end
           
           if ~fleq(G.ioDelay,0)
              G = ss(G); 
           end
           
           % Initial controller
           cin    = Z.initialParams;
           Gcinit = oustapid(cin(1),cin(2),cin(4),cin(3),cin(5),w(1),w(2),N,type);
           
           % Final controller
           Gcfinal = oustapid(Kp,Ki,lam,Kd,mu,w(1),w(2),N,type);
           
           ig1 = []; ig2 = [];
           
           if isempty(usesim)
               
               % Use lsim() for system simulation
  
               % Check systems to be proper
               if get(handles.chkZeroCancelation, 'Value')
                   G = toproper(G, w(2));
                   Gcinit= toproper(Gcinit, w(2));
                   Gcfinal = toproper(Gcfinal, w(2));
               end
               
               % Get initial step response
               y_init = lsim(feedback(G*Gcinit,1),sp*ones(length(t),1),t);

               % Get new step response
               y = lsim(feedback(G*Gcfinal,1),sp*ones(length(t),1),t);
               
               % Set time vectors for compare
               t1 = t;
               t2 = t;
               
               ref = sp;
               
           else
              
               % Use Simulink for simulation
               sopt = fpid_simopt(usesim, maxtime, dtmin, dtmax, ulim, ...
                                 sp, [w(1) w(2)], N, type, cancelzero);
               
               % Initial system
               [y_init, ig1, t1, ref] = fpid_optimize_sim([cin(1),cin(2),cin(4),cin(3),cin(5)], G, sopt);
               
               % Final system
               [y, ig2, t2, ref] = fpid_optimize_sim([Kp,Ki,lam,Kd,mu], G, sopt);
               
           end
           
           legendTimeDomainYtText = {'Initial response'};
           legendTimeDomainUtText = {'Initial control law'};
           legendFreqDomainText   = {'Initial frequency response'};
           
           % Legend for the signals
           if ~get(handles.chkSimulateOnly, 'Value')
               legendTimeDomainYtText{end+1} = 'Post-optimization response';
               legendTimeDomainUtText{end+1} = 'Post-optimization control law';
               legendFreqDomainText{end+1}   = 'Post-optimization frequency response';
           else
               legendTimeDomainYtText = strrep(legendTimeDomainYtText, 'Initial', 'Current');
               legendTimeDomainUtText = strrep(legendTimeDomainUtText, 'Initial', 'Current');
               legendFreqDomainText   = strrep(legendFreqDomainText, 'Initial', 'Current');
           end
           
           % Closed-loop simulation: Plot control law if data present
           h0 = figure();
           hold on;
           set(h0, 'NumberTitle', 'off');
           set(h0, 'Name', 'Controller optimization results');
           
           if (isempty(ig1) || isempty(ig2))
               h01 = plot(t1, y_init, 'b-', 'Linewidth', 2);
               h02 = plot(t2, y, 'g--', 'Linewidth', 2);
               h03 = plot(t2, ref, 'r:');
               legend(legendTimeDomainYtText, 'Location', 'Best');
               grid;
           else
               % System response
               subplot(2,1,1);
               hold on;
               h01 = plot(t1, y_init, 'b-', 'Linewidth', 2);
               h02 = plot(t2, y, 'g--', 'Linewidth', 2);
               h03 = plot(t2, ref, 'r:');
               the_clearance = 0.05*abs(max(y_init)-min(y_init));
               y_lim = get(gca, 'ylim');
               ylim([y_lim(1)-the_clearance y_lim(2)+the_clearance]);
               legend(legendTimeDomainYtText, 'Location', 'Best');
               ylabel('Amplitude');
               grid;
               
               % Control law
               subplot(2,1,2);
               hold on;
               h011 = plot(t1, ig1, 'b-', 'Linewidth', 2);
               h012 = plot(t2, ig2, 'g--', 'Linewidth', 2);
               legend(legendTimeDomainUtText, 'Location', 'Best');
               ylabel('Control law u(t)');
               xlabel('Time [s]');
               the_clearance = 0.05*abs(max(ig2)-min(ig2));
               y_lim = get(gca, 'ylim');
               ylim([y_lim(1)-the_clearance y_lim(2)+the_clearance]);
               grid;
           end
           
           if get(handles.chkSimulateOnly, 'Value')
               delete(findall(h0, 'LineStyle', '--'));
           end
           
           % Frequency range for all frequency response plots
           w_range = logspace(get10exp(w(1)),get10exp(w(2)),1000);
           
           if get(handles.chkGMEnable, 'Value') || ...
              get(handles.chkEnableRobustness, 'Value')
               % Open-loop Bode diagram
               h1  = figure();
               hold on;
               h11 = bodeplot(G*Gcinit, 'b-', w_range);
               h12 = bodeplot(G*Gcfinal, 'g--', w_range);
               
               set(h1, 'NumberTitle', 'off');
               set(h1, 'Name', 'Open-loop frequency response comparison');
               set(findall(gcf,'type','line'),'linewidth',2);
               legend(legendFreqDomainText, 'Location', 'Best');
               grid;
               
               if get(handles.chkSimulateOnly, 'Value')
                   set(findall(h1, 'LineStyle', '--'), 'Visible', 'Off');
               end
           end
           
           if get(handles.chkSensitivity, 'Value')
               % Complementary sensitivity function
               h2  = figure();
               hold on;
               h21 = bodeplot(G*Gcinit / (1+G*Gcinit), 'b-', w_range);
               h22 = bodeplot(G*Gcfinal / (1+G*Gcfinal), 'g--', w_range);
               
               set(h2, 'NumberTitle', 'off');
               set(h2, 'Name', 'Complementary sensitivity function T(jw) comparison');
               set(findall(gcf,'type','line'),'linewidth',2);
               legend(legendFreqDomainText, 'Location', 'Best');
               grid;
               
               if get(handles.chkSimulateOnly, 'Value')
                   set(findall(h2, 'LineStyle', '--'), 'Visible', 'Off');
               end
               
               % Sensitivity function
               h3  = figure();
               hold on;
               h31 = bodeplot(1/(1+G*Gcinit), 'b-', w_range);
               h32 = bodeplot(1/(1+G*Gcfinal), 'g--', w_range);
               set(h3, 'NumberTitle', 'off');
               set(h3, 'Name', 'Sensitivity function S(jw) comparison');
               set(findall(gcf,'type','line'),'linewidth',2);
               legend(legendFreqDomainText, 'Location', 'Best');
               grid;
               
               if get(handles.chkSimulateOnly, 'Value')
                   set(findall(h3, 'LineStyle', '--'), 'Visible', 'Off');
               end
               
           end

        end
 
    else
        
        errordlg('System no longer in workspace or invalid!');
        
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
    set(ft.txtLambda,'String',get(handles.txtLam,'String'));
    set(ft.txtMu,'String',get(handles.txtMu,'String'));

% --- Executes on button press in chkSimulate.
function chkSimulate_Callback(hObject, eventdata, handles)
% hObject    handle to chkSimulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkSimulate



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



function txtAppOrder_Callback(hObject, eventdata, handles)
% hObject    handle to txtAppOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAppOrder as text
%        str2double(get(hObject,'String')) returns contents of txtAppOrder as a double


% --- Executes during object creation, after setting all properties.
function txtAppOrder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAppOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtKpMin_Callback(hObject, eventdata, handles)
% hObject    handle to txtKpMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtKpMin as text
%        str2double(get(hObject,'String')) returns contents of txtKpMin as a double


% --- Executes during object creation, after setting all properties.
function txtKpMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtKpMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtKpMax_Callback(hObject, eventdata, handles)
% hObject    handle to txtKpMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtKpMax as text
%        str2double(get(hObject,'String')) returns contents of txtKpMax as a double


% --- Executes during object creation, after setting all properties.
function txtKpMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtKpMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtKiMin_Callback(hObject, eventdata, handles)
% hObject    handle to txtKiMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtKiMin as text
%        str2double(get(hObject,'String')) returns contents of txtKiMin as a double


% --- Executes during object creation, after setting all properties.
function txtKiMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtKiMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtKiMax_Callback(hObject, eventdata, handles)
% hObject    handle to txtKiMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtKiMax as text
%        str2double(get(hObject,'String')) returns contents of txtKiMax as a double


% --- Executes during object creation, after setting all properties.
function txtKiMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtKiMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtLamMin_Callback(hObject, eventdata, handles)
% hObject    handle to txtLamMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtLamMin as text
%        str2double(get(hObject,'String')) returns contents of txtLamMin as a double


% --- Executes during object creation, after setting all properties.
function txtLamMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtLamMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtLamMax_Callback(hObject, eventdata, handles)
% hObject    handle to txtLamMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtLamMax as text
%        str2double(get(hObject,'String')) returns contents of txtLamMax as a double


% --- Executes during object creation, after setting all properties.
function txtLamMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtLamMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtKdMin_Callback(hObject, eventdata, handles)
% hObject    handle to txtKdMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtKdMin as text
%        str2double(get(hObject,'String')) returns contents of txtKdMin as a double


% --- Executes during object creation, after setting all properties.
function txtKdMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtKdMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtKdMax_Callback(hObject, eventdata, handles)
% hObject    handle to txtKdMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtKdMax as text
%        str2double(get(hObject,'String')) returns contents of txtKdMax as a double


% --- Executes during object creation, after setting all properties.
function txtKdMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtKdMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txtMuMin_Callback(hObject, eventdata, handles)
% hObject    handle to txtMuMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtMuMin as text
%        str2double(get(hObject,'String')) returns contents of txtMuMin as a double


% --- Executes during object creation, after setting all properties.
function txtMuMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtMuMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtMuMax_Callback(hObject, eventdata, handles)
% hObject    handle to txtMuMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtMuMax as text
%        str2double(get(hObject,'String')) returns contents of txtMuMax as a double


% --- Executes during object creation, after setting all properties.
function txtMuMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtMuMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtGm_Callback(hObject, eventdata, handles)
% hObject    handle to txtGm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGm as text
%        str2double(get(hObject,'String')) returns contents of txtGm as a double


% --- Executes during object creation, after setting all properties.
function txtGm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtPm_Callback(hObject, eventdata, handles)
% hObject    handle to txtPm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtPm as text
%        str2double(get(hObject,'String')) returns contents of txtPm as a double


% --- Executes during object creation, after setting all properties.
function txtPm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtPm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkStrict.
function chkStrict_Callback(hObject, eventdata, handles)
% hObject    handle to chkStrict (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkStrict


% --- Executes on button press in chkLimitIter.
function chkLimitIter_Callback(hObject, eventdata, handles)
% hObject    handle to chkLimitIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkLimitIter
    
    maxIterState = {'off', 'on'};
    set(handles.txtLimitIter,'Enable',maxIterState{get(hObject,'Value')+1});


function txtLimitIter_Callback(hObject, eventdata, handles)
% hObject    handle to txtLimitIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtLimitIter as text
%        str2double(get(hObject,'String')) returns contents of txtLimitIter as a double


% --- Executes during object creation, after setting all properties.
function txtLimitIter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtLimitIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menuFix.
function menuFix_Callback(hObject, eventdata, handles)
% hObject    handle to menuFix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menuFix contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuFix
    
    % Toggle parameter input on/off
    switch get(hObject,'Value')
        case 1
            expEnable = 'on';
            gainEnable = 'on';
        case 2
            expEnable = 'off';
            gainEnable = 'on';
        case 3
            expEnable = 'on';
            gainEnable = 'off';
    end
    
    % Exponents
    set(handles.txtLam,    'Enable', expEnable);
    set(handles.txtMu,     'Enable', expEnable);
    set(handles.txtLamMin, 'Enable', expEnable);
    set(handles.txtMuMin,  'Enable', expEnable);
    set(handles.txtLamMax, 'Enable', expEnable);
    set(handles.txtMuMax,  'Enable', expEnable);
    
    % Gains
    set(handles.txtKp,    'Enable', gainEnable);
    set(handles.txtKpMin, 'Enable', gainEnable);
    set(handles.txtKpMax, 'Enable', gainEnable);
    
	set(handles.txtKi,    'Enable', gainEnable);
    set(handles.txtKiMin, 'Enable', gainEnable);
    set(handles.txtKiMax, 'Enable', gainEnable);
    
	set(handles.txtKd,    'Enable', gainEnable);
    set(handles.txtKdMin, 'Enable', gainEnable);
    set(handles.txtKdMax, 'Enable', gainEnable);
    

% --- Executes during object creation, after setting all properties.
function menuFix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuFix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtPlantType_Callback(hObject, eventdata, handles)
% hObject    handle to txtPlantType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtPlantType as text
%        str2double(get(hObject,'String')) returns contents of txtPlantType as a double


% --- Executes during object creation, after setting all properties.
function txtPlantType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtPlantType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkPlantType(handles)
    
    plantName = get(handles.txtPlant, 'String');
    
    if ~isempty(plantName)
        
        sys = evalin('base', ['class(' plantName ')']);

        switch sys
            case 'fotf'
            case 'tf'
            case 'zpk'
            case 'ss'
            otherwise
                sys = '';
        end
        
    else
        
        sys = '';
        
    end
    
    set(handles.txtPlantType, 'String', sys);



function txtSP_Callback(hObject, eventdata, handles)
% hObject    handle to txtSP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtSP as text
%        str2double(get(hObject,'String')) returns contents of txtSP as a double


% --- Executes during object creation, after setting all properties.
function txtSP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkControlLimits.
function chkControlLimits_Callback(hObject, eventdata, handles)
% hObject    handle to chkControlLimits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkControlLimits
if get(hObject, 'Value')
   % Turn on gain & phase margin specifications
   set(handles.txtUMin, 'Enable', 'On');
   set(handles.txtUMax, 'Enable', 'On');
else
   % Else turn them off
   set(handles.txtUMin, 'Enable', 'Off');
   set(handles.txtUMax, 'Enable', 'Off');
end


function txtUMin_Callback(hObject, eventdata, handles)
% hObject    handle to txtUMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtUMin as text
%        str2double(get(hObject,'String')) returns contents of txtUMin as a double


% --- Executes during object creation, after setting all properties.
function txtUMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtUMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtUMax_Callback(hObject, eventdata, handles)
% hObject    handle to txtUMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtUMax as text
%        str2double(get(hObject,'String')) returns contents of txtUMax as a double


% --- Executes during object creation, after setting all properties.
function txtUMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtUMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkSensitivity.
function chkSensitivity_Callback(hObject, eventdata, handles)
% hObject    handle to chkSensitivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkSensitivity
if get(hObject, 'Value')
   % Turn on gain & phase margin specifications
   set(handles.txtTdb, 'Enable', 'On');
   set(handles.txtSdb, 'Enable', 'On');
   set(handles.txtWt, 'Enable', 'On');
   set(handles.txtWs, 'Enable', 'On');
else
   % Else turn them off
   set(handles.txtTdb, 'Enable', 'Off');
   set(handles.txtSdb, 'Enable', 'Off');
   set(handles.txtWt, 'Enable', 'Off');
   set(handles.txtWs, 'Enable', 'Off');
end


function txtTdb_Callback(hObject, eventdata, handles)
% hObject    handle to txtTdb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTdb as text
%        str2double(get(hObject,'String')) returns contents of txtTdb as a double


% --- Executes during object creation, after setting all properties.
function txtTdb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTdb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtWt_Callback(hObject, eventdata, handles)
% hObject    handle to txtWt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtWt as text
%        str2double(get(hObject,'String')) returns contents of txtWt as a double


% --- Executes during object creation, after setting all properties.
function txtWt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtWt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtSdb_Callback(hObject, eventdata, handles)
% hObject    handle to txtSdb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtSdb as text
%        str2double(get(hObject,'String')) returns contents of txtSdb as a double


% --- Executes during object creation, after setting all properties.
function txtSdb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSdb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtWs_Callback(hObject, eventdata, handles)
% hObject    handle to txtWs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtWs as text
%        str2double(get(hObject,'String')) returns contents of txtWs as a double


% --- Executes during object creation, after setting all properties.
function txtWs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtWs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkGMEnable.
function chkGMEnable_Callback(hObject, eventdata, handles)
% hObject    handle to chkGMEnable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject, 'Value')
   % Turn on gain & phase margin specifications
   set(handles.txtGm, 'Enable', 'On');
   set(handles.chkGmExact, 'Enable', 'On');
   set(handles.txtPm, 'Enable', 'On');
   set(handles.chkPmExact, 'Enable', 'On');
else
   % Else turn them off
   set(handles.txtGm, 'Enable', 'Off');
   set(handles.chkGmExact, 'Enable', 'Off');
   set(handles.txtPm, 'Enable', 'Off');
   set(handles.chkPmExact, 'Enable', 'Off');
end

% Hint: get(hObject,'Value') returns toggle state of chkGMEnable


% --- Executes on button press in chkZeroCancelation.
function chkZeroCancelation_Callback(hObject, eventdata, handles)
% hObject    handle to chkZeroCancelation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkZeroCancelation



function txtSimResolution_Callback(hObject, eventdata, handles)
% hObject    handle to txtSimResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtSimResolution as text
%        str2double(get(hObject,'String')) returns contents of txtSimResolution as a double


% --- Executes during object creation, after setting all properties.
function txtSimResolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSimResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnReset.
function btnReset_Callback(hObject, eventdata, handles)
% hObject    handle to btnReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    set_gains = false;
    set_exponents = false;
    
    reset_to = get(handles.txtReset, 'String');
    
    switch(get(handles.menuReset, 'Value'))
        case 1
            set_gains = true;
        case 2
            set_exponents = true;
        case 3
            set_gains = true;
            set_exponents = true;
    end
    
    if set_gains
        set(handles.txtKp, 'String', reset_to);
        set(handles.txtKi, 'String', reset_to);
        set(handles.txtKd, 'String', reset_to);
    end
    
    if set_exponents
        set(handles.txtLam, 'String', reset_to);
        set(handles.txtMu, 'String', reset_to);
    end


function txtReset_Callback(hObject, eventdata, handles)
% hObject    handle to txtReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtReset as text
%        str2double(get(hObject,'String')) returns contents of txtReset as a double


% --- Executes during object creation, after setting all properties.
function txtReset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menuReset.
function menuReset_Callback(hObject, eventdata, handles)
% hObject    handle to menuReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menuReset contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuReset


% --- Executes during object creation, after setting all properties.
function menuReset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menuOptimizer.
function menuOptimizer_Callback(hObject, eventdata, handles)
% hObject    handle to menuOptimizer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menuOptimizer contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuOptimizer



% --- Executes during object creation, after setting all properties.
function menuOptimizer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuOptimizer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menuOptimizerAlgorithm.
function menuOptimizerAlgorithm_Callback(hObject, eventdata, handles)
% hObject    handle to menuOptimizerAlgorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Remove strict option if anything else than 'nelder-mead' is selected
    switch(get(hObject, 'Value'));
        case 1
            set(handles.chkStrict, 'Enable', 'On');
        case {2, 3, 4}
            set(handles.chkStrict, 'Enable', 'Off');
    end
    

% Hints: contents = cellstr(get(hObject,'String')) returns menuOptimizerAlgorithm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuOptimizerAlgorithm


% --- Executes during object creation, after setting all properties.
function menuOptimizerAlgorithm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuOptimizerAlgorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txtMaxTime_Callback(hObject, eventdata, handles)
% hObject    handle to txtMaxTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtMaxTime as text
%        str2double(get(hObject,'String')) returns contents of txtMaxTime as a double


% --- Executes during object creation, after setting all properties.
function txtMaxTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtMaxTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtDtMin_Callback(hObject, eventdata, handles)
% hObject    handle to txtDtMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtDtMin as text
%        str2double(get(hObject,'String')) returns contents of txtDtMin as a double


% --- Executes during object creation, after setting all properties.
function txtDtMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtDtMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtDtMax_Callback(hObject, eventdata, handles)
% hObject    handle to txtDtMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtDtMax as text
%        str2double(get(hObject,'String')) returns contents of txtDtMax as a double


% --- Executes during object creation, after setting all properties.
function txtDtMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtDtMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkEnableRobustness.
function chkEnableRobustness_Callback(hObject, eventdata, handles)
% hObject    handle to chkEnableRobustness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkEnableRobustness
if get(hObject, 'Value')
   % Turn on gain & phase margin specifications
   set(handles.txtW, 'Enable', 'On');
   set(handles.chkWh, 'Enable', 'On');
   if get(handles.chkWh, 'Value')
        set(handles.txtRange, 'Enable', 'On');
   end
else
   % Else turn them off
   set(handles.chkWh, 'Enable', 'Off');
   set(handles.txtW, 'Enable', 'Off');
   set(handles.txtRange, 'Enable', 'Off');
end


function txtW_Callback(hObject, eventdata, handles)
% hObject    handle to txtW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtW as text
%        str2double(get(hObject,'String')) returns contents of txtW as a double


% --- Executes during object creation, after setting all properties.
function txtW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtRange_Callback(hObject, eventdata, handles)
% hObject    handle to txtRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtRange as text
%        str2double(get(hObject,'String')) returns contents of txtRange as a double


% --- Executes during object creation, after setting all properties.
function txtRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkUseSimulink.
function chkUseSimulink_Callback(hObject, eventdata, handles)
% hObject    handle to chkUseSimulink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkUseSimulink
    
    % Saturation specifications only available if
    % Simulink-based simulation is used
    if get(hObject, 'Value')
        set(handles.chkControlLimits, 'Enable', 'On');
        set(handles.chkWeight, 'Enable', 'On');
        if (get(handles.chkControlLimits, 'Value'))
            set(handles.txtUMin, 'Enable', 'On');
            set(handles.txtUMax, 'Enable', 'On');
        end
        if (get(handles.chkWeight, 'Value'))
            set(handles.txtWeight, 'Enable', 'On');
        end
        set(handles.menuSimulinkModel, 'Enable', 'On');
        if get(handles.menuSimulinkModel, 'Value') > 1
            set(handles.btnSimEdit, 'Enable', 'On');
        end
        set(handles.btnSimNew, 'Enable', 'On');
    else
        set(handles.chkControlLimits, 'Enable', 'Off');
        set(handles.txtUMin, 'Enable', 'Off');
        set(handles.txtUMax, 'Enable', 'Off');
        set(handles.chkWeight, 'Enable', 'Off');
        set(handles.txtWeight, 'Enable', 'Off');
        set(handles.menuSimulinkModel, 'Enable', 'Off');
        set(handles.btnSimEdit, 'Enable', 'Off');
        set(handles.btnSimNew, 'Enable', 'Off');
    end
    
    % Refresh model list
    refreshSimulinkModels(handles);
    
    


% --- Executes on button press in chkDisableWarnings.
function chkDisableWarnings_Callback(hObject, eventdata, handles)
% hObject    handle to chkDisableWarnings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkDisableWarnings

% Disable warnings, if desired
if get(handles.chkDisableWarnings, 'Value')
    warns = warning('query', 'all');
    opt = get(handles.output, 'UserData');
    opt.warnings_state = warns;
    set(handles.output, 'UserData', opt);
    warning off all;
else
    opt = get(handles.output, 'UserData');
    [ex, warns] = cfieldexists(opt, 'warnings_state'); %#ok<ASGLU>
    if ~isempty(warns), warning(warns); end
end

% --- Executes on button press in chkPmExact.
function chkPmExact_Callback(hObject, eventdata, handles)
% hObject    handle to chkPmExact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkPmExact


% --- Executes on button press in chkGmExact.
function chkGmExact_Callback(hObject, eventdata, handles)
% hObject    handle to chkGmExact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkGmExact


% --- Executes on selection change in menuSimulinkModel.
function menuSimulinkModel_Callback(hObject, eventdata, handles)
% hObject    handle to menuSimulinkModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menuSimulinkModel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuSimulinkModel
if get(handles.menuSimulinkModel, 'Value') > 1
    set(handles.btnSimEdit, 'Enable', 'On');
else
    set(handles.btnSimEdit, 'Enable', 'Off');
end


% --- Executes during object creation, after setting all properties.
function menuSimulinkModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuSimulinkModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnSimNew.
function btnSimNew_Callback(hObject, eventdata, handles)
% hObject    handle to btnSimNew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
name        = 'Default system copy creation';
desc        = 'The default system will be copied into current working directory under name';
prompt      = {'Model name:'};
defaultName = {'fpid_optimize_'};

options.Resize      = 'on';
options.WindowStyle = 'normal';

crNew = inputdlg(prompt,name,1,defaultName,options);

if ~isempty(crNew)
    
    % Get file name
    newFileName = crNew{1};
    
    try
        % Check extension
        [ig1, ig2, ext] = fileparts(newFileName);
        if isempty(ext), newFileName = [newFileName '.mdl']; end
    
        % Copy default model to new file
        copyfile(which('fpid_optimize_model.mdl'), [pwd filesep newFileName], 'f');
    
    catch e
       
        % Failed, do nothing
        errordlg(e.message);
       
    end
    
    % Refresh model list
    refreshSimulinkModels(handles);
    
    % Find the newly created model
    allModels = get(handles.menuSimulinkModel, 'String');
    for i=1:length(allModels)
       if strcmpi(allModels{i},newFileName)
           set(handles.menuSimulinkModel, 'Value', i);
           set(handles.btnSimEdit, 'Enable', 'On');
       end
    end
end
    

function refreshSimulinkModels(handles)
    % Update the available model list
    allCwdFiles = what();           % Get models from current
    allModels   = allCwdFiles.mdl;  % working directory                                
    modelList = ['default'; allModels(:)];
    
    % Check for index "overflow"
    if numel(modelList) < get(handles.menuSimulinkModel, 'Value')
        set(handles.menuSimulinkModel, 'Value', 1);
    end
    
    set(handles.menuSimulinkModel, 'String', modelList);


% --- Executes on button press in btnSimEdit.
function btnSimEdit_Callback(hObject, eventdata, handles)
% hObject    handle to btnSimEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    toEdit = get(handles.menuSimulinkModel, 'String');
    open_system(which(toEdit{get(handles.menuSimulinkModel, 'Value')}));
catch e
    errordlg(e.message);
    refreshSimulinkModels(handles);
end


% --- Executes on button press in chkWeight.
function chkWeight_Callback(hObject, eventdata, handles)
% hObject    handle to chkWeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkWeight
if get(handles.chkWeight, 'Value')
    set(handles.txtWeight, 'Enable', 'On');
else
    set(handles.txtWeight, 'Enable', 'Off');
end


function txtWeight_Callback(hObject, eventdata, handles)
% hObject    handle to txtWeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtWeight as text
%        str2double(get(hObject,'String')) returns contents of txtWeight as a double


% --- Executes during object creation, after setting all properties.
function txtWeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtWeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkWh.
function chkWh_Callback(hObject, eventdata, handles)
% hObject    handle to chkWh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkWh
if get(hObject,'Value') && get(handles.chkEnableRobustness,'Value')
    set(handles.txtRange,'Enable','On');
else
    set(handles.txtRange,'Enable','Off');
end


% --------------------------------------------------------------------
function menuTools_Callback(hObject, eventdata, handles)
% hObject    handle to menuTools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuOptions_Callback(hObject, eventdata, handles)
% hObject    handle to menuOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
config = fomcon('config');
fpid_config = struct;
fpid_config.FPID_Optimizer = config.FPID_Optimizer;
[h, new_fpid_config] = propertiesGUI([], fpid_config);

% Save to workspace
if ~isempty(new_fpid_config)
    config.FPID_Optimizer = new_fpid_config.FPID_Optimizer;
    fomcon('config', config);
end


% --------------------------------------------------------------------
function menuFile_Callback(hObject, eventdata, handles)
% hObject    handle to menuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuLoadConfig_Callback(hObject, eventdata, handles)
% hObject    handle to menuLoadConfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,path] = uigetfile();
if ~isequal(filename, 0) && ~isequal(path, 0)
    config_struct   = load([path filename]);
    config = config_struct.FPID_Optimizer_GUI_config;
    
    % Handles shortcut
    h = handles;
    
    % Fetch the main structures
    PlantModel = config.PlantModel;
    FPIDParams = config.FPIDParams;
    Simulation = config.Simulation;
    Optimization = config.Optimization;
    
    % Set all the parameters
    
    % Plant model
    set(h.txtPlant, 'String', PlantModel.LTIsystem);
    set(h.txtPlantType, 'String', PlantModel.Type);
    set(h.menuApprox, 'Value', PlantModel.ApproxType);
    set(h.txtFreqRange, 'String', PlantModel.ApproxFreqRange);
    set(h.txtAppOrder, 'String', PlantModel.ApproxOrder);
    set(h.chkZeroCancelation, 'Value', PlantModel.ZeroCancelation);
    
    % Fractional PID parameters
    set(h.menuFix, 'Value', FPIDParams.ParameterFix); % Callback
    menuFix_Callback(h.menuFix,[],h);
    set(h.txtKp, 'String', FPIDParams.Kp);
    set(h.txtKpMin, 'String', FPIDParams.KpMin);
    set(h.txtKpMax, 'String', FPIDParams.KpMax);
    set(h.txtKi, 'String', FPIDParams.Ki);
    set(h.txtKiMin, 'String', FPIDParams.KiMin);
    set(h.txtKiMax, 'String', FPIDParams.KiMax);
    set(h.txtKd, 'String', FPIDParams.Kd);
    set(h.txtKdMin, 'String', FPIDParams.KdMin);
    set(h.txtKdMax, 'String', FPIDParams.KdMax);
    set(h.menuReset, 'Value', FPIDParams.Reset);
    set(h.txtReset, 'String', FPIDParams.ResetVal);
    
    set(h.txtLam, 'String', FPIDParams.Lam);
    set(h.txtLamMin, 'String', FPIDParams.LamMin);
    set(h.txtLamMax, 'String', FPIDParams.LamMax);
    
    set(h.txtMu, 'String', FPIDParams.Mu);
    set(h.txtMuMin, 'String', FPIDParams.MuMin);
    set(h.txtMuMax, 'String', FPIDParams.MuMax);
    
    % Simulation parameters
    set(h.txtMaxTime, 'String', Simulation.MaxTime);
    set(h.txtDtMin, 'String', Simulation.DtMin);
    set(h.txtDtMax, 'String', Simulation.DtMax);
    set(h.chkUseSimulink, 'Value', Simulation.Simulink); % Callback
    chkUseSimulink_Callback(h.chkUseSimulink,[],h);
    set(h.chkDisableWarnings, 'Value', Simulation.DisableWarnings);
    if ~Simulation.UseStandardModel
        warndlg('Do not forget to select the correct Simulink model!', ...
            'Non-default model warning', ...
            'modal');
    end
    
    % Optimization options
    set(h.menuOptimizerAlgorithm, 'Value', Optimization.Algorithm); % Callback
    menuOptimizerAlgorithm_Callback(h.menuOptimizerAlgorithm,[],h);
    set(h.menuMetric, 'Value', Optimization.Metric);
    
    set(h.chkGMEnable, 'Value', Optimization.GmPm.Enable); % Callback
    chkGMEnable_Callback(h.chkGMEnable, [], h);
    set(h.txtGm, 'String', Optimization.GmPm.Gm);
    set(h.chkGmExact, 'Value', Optimization.GmPm.GmExact);
    set(h.txtPm, 'String', Optimization.GmPm.Pm);
    set(h.chkPmExact, 'Value', Optimization.GmPm.PmExact);
    
    set(h.chkSensitivity, 'Value', Optimization.SensFun.Enable); % Callback
    chkSensitivity_Callback(h.chkSensitivity, [], h);
    set(h.txtTdb, 'String', Optimization.SensFun.T);
    set(h.txtWt, 'String', Optimization.SensFun.Wt);
    set(h.txtSdb, 'String', Optimization.SensFun.S);
    set(h.txtWs, 'String', Optimization.SensFun.Ws);
    
    set(h.chkEnableRobustness, 'Value', Optimization.GainVar.Enable); % Callback
    chkEnableRobustness_Callback(h.chkEnableRobustness, [], h);
    set(h.txtW, 'String', Optimization.GainVar.Wc);
    set(h.chkWh, 'Value', Optimization.GainVar.WhEnable); % Callback
    chkWh_Callback(h.chkWh, [], h);
    set(h.txtRange, 'String', Optimization.GainVar.Wh);
    
    set(h.chkControlLimits, 'Value', Optimization.ControlLaw.Enable); % Callback
    chkControlLimits_Callback(h.chkControlLimits, [], h);
    set(h.chkWeight, 'Value', Optimization.ControlLaw.EnableWeight); % Callback
    chkWeight_Callback(h.chkWeight, [], h);
    set(h.txtWeight, 'String', Optimization.ControlLaw.Weight);
    set(h.txtUMin, 'String', Optimization.ControlLaw.UMin);
    set(h.txtUMax, 'String', Optimization.ControlLaw.UMax);
    
    set(h.txtSP, 'String', Optimization.SetPoint);
    set(h.chkStrict, 'Value', Optimization.ForceStrict);
    set(h.chkSimulate, 'Value', Optimization.GenerateReport);
    set(h.chkLimitIter, 'Value', Optimization.EnableIterationLimit); % Callback
    chkLimitIter_Callback(h.chkLimitIter, [], h);
    set(h.txtLimitIter, 'String', Optimization.IterationLimit);
    
    % New options after initial release. Must provide
    % a backwards-compatibility mechanism
    
    % FOMCON global FPID optimizer options
    conf = fomcon('config'); % Make sure the configuration is in the base workspace
    
    [fe, f] = cfieldexists(config, 'FOMCON.Time_domain_computations.Performance_index_weight'); %#ok<ASGLU>
    if ~isempty(f)
        conf.FPID_Optimizer.Time_domain_computations.Performance_index_weight = f;
    else
        warning('FPIDOPTIM:ConfigurationParameterMissing', ...
                'Parameter ''Performance_index_weight'' is missing and was not assigned');
    end
    
    [fe, f] = cfieldexists(config, 'FOMCON.Frequency_domain_computations.Num_points_phase'); %#ok<ASGLU>
    if ~isempty(f)
        conf.FPID_Optimizer.Frequency_domain_computations.Num_points_phase = f;
    else
        warning('FPIDOPTIM:ConfigurationParameterMissing', ...
                'Parameter ''Num_points_phase'' is missing and was not assigned');
    end
    
    [fe, f] = cfieldexists(config, 'FOMCON.Frequency_domain_computations.Phase_comp_weight'); %#ok<ASGLU>
    if ~isempty(f)
        conf.FPID_Optimizer.Frequency_domain_computations.Phase_comp_weight = f;
    else
        warning('FPIDOPTIM:ConfigurationParameterMissing', ...
                'Parameter ''Phase_comp_weight'' is missing and was not assigned');
    end
    
    [fe, f] = cfieldexists(config, 'FOMCON.Frequency_domain_computations.W_cg_comp_weight'); %#ok<ASGLU>
    if ~isempty(f)
        conf.FPID_Optimizer.Frequency_domain_computations.W_cg_comp_weight = f;
    else
        warning('FPIDOPTIM:ConfigurationParameterMissing', ...
                'Parameter ''W_cg_comp_weight'' is missing and was not assigned');
    end
    
    % Save the FOMCON structure
    fomcon('config', conf);
    
    % Do the "Simulate only" callback
    chkSimulateOnly_Callback(h.chkSimulateOnly, [], h);
    
    % Do the warnings callback
    chkDisableWarnings_Callback(h.chkDisableWarnings, [], h);
end

% --------------------------------------------------------------------
function menuSaveConfig_Callback(hObject, eventdata, handles)
% hObject    handle to menuSaveConfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Open file for writing
[filename,path] = uiputfile('*.mat');

if ~isequal(filename, 0) && ~isequal(path, 0)

    % Check file extension, add if necessary
    if ~strcmpi(filename(end-2:end),'mat')
        filename = [filename '.mat'];
    end

    % Create the file
    config_struct = matfile([path filename], 'Writable', true);

    % Handles shortcut
    h = handles;

    % Build the configuration structure

    % Plant model
    PlantModel = struct;
    PlantModel.LTIsystem = get(h.txtPlant, 'String');
    PlantModel.Type = get(h.txtPlantType, 'String');
    PlantModel.ApproxType = get(h.menuApprox, 'Value');
    PlantModel.ApproxFreqRange = get(h.txtFreqRange, 'String');
    PlantModel.ApproxOrder = get(h.txtAppOrder, 'String');
    PlantModel.ZeroCancelation = get(h.chkZeroCancelation, 'Value');

    % Fractional PID parameters
    FPIDParams = struct;
    FPIDParams.ParameterFix = get(h.menuFix, 'Value'); % Callback
    FPIDParams.Kp = get(h.txtKp, 'String');
    FPIDParams.KpMin = get(h.txtKpMin, 'String');
    FPIDParams.KpMax = get(h.txtKpMax, 'String');
    FPIDParams.Ki = get(h.txtKi, 'String');
    FPIDParams.KiMin = get(h.txtKiMin, 'String');
    FPIDParams.KiMax = get(h.txtKiMax, 'String');
    FPIDParams.Kd = get(h.txtKd, 'String');
    FPIDParams.KdMin = get(h.txtKdMin, 'String');
    FPIDParams.KdMax = get(h.txtKdMax, 'String');
    
    FPIDParams.Lam = get(h.txtLam, 'String');
    FPIDParams.LamMin = get(h.txtLamMin, 'String');
    FPIDParams.LamMax = get(h.txtLamMax, 'String');
    
    FPIDParams.Mu = get(h.txtMu, 'String');
    FPIDParams.MuMin = get(h.txtMuMin, 'String');
    FPIDParams.MuMax = get(h.txtMuMax, 'String');
    
    FPIDParams.Reset = get(h.menuReset, 'Value');
    FPIDParams.ResetVal = get(h.txtReset, 'String');

    % Simulation parameters
    Simulation = struct;
    Simulation.MaxTime = get(h.txtMaxTime, 'String');
    Simulation.DtMin = get(h.txtDtMin, 'String');
    Simulation.DtMax = get(h.txtDtMax, 'String');
    Simulation.Simulink = get(h.chkUseSimulink, 'Value'); % Callback
    Simulation.DisableWarnings = get(h.chkDisableWarnings, 'Value');
    if get(h.menuSimulinkModel, 'Value')>1
        Simulation.UseStandardModel = false;
    else
        Simulation.UseStandardModel = true;
    end

    % Optimization options
    Optimization = struct;
    Optimization.Algorithm = get(h.menuOptimizerAlgorithm, 'Value'); % Callback
    Optimization.Metric = get(h.menuMetric, 'Value');

    Optimization.GmPm.Enable = get(h.chkGMEnable, 'Value'); % Callback
    Optimization.GmPm.Gm = get(h.txtGm, 'String');
    Optimization.GmPm.GmExact = get(h.chkGmExact, 'Value');
    Optimization.GmPm.Pm = get(h.txtPm, 'String');
    Optimization.GmPm.PmExact = get(h.chkPmExact, 'Value');

    Optimization.SensFun.Enable = get(h.chkSensitivity, 'Value'); % Callback
    Optimization.SensFun.T = get(h.txtTdb, 'String');
    Optimization.SensFun.Wt = get(h.txtWt, 'String');
    Optimization.SensFun.S = get(h.txtSdb, 'String');
    Optimization.SensFun.Ws = get(h.txtWs, 'String');

    Optimization.GainVar.Enable = get(h.chkEnableRobustness, 'Value'); % Callback
    Optimization.GainVar.Wc = get(h.txtW, 'String');
    Optimization.GainVar.WhEnable = get(h.chkWh, 'Value');
    Optimization.GainVar.Wh = get(h.txtRange, 'String');

    Optimization.ControlLaw.Enable = get(h.chkControlLimits, 'Value'); % Callback
    Optimization.ControlLaw.EnableWeight = get(h.chkWeight, 'Value'); % Callback
    Optimization.ControlLaw.Weight = get(h.txtWeight, 'String');
    Optimization.ControlLaw.UMin = get(h.txtUMin, 'String');
    Optimization.ControlLaw.UMax = get(h.txtUMax, 'String');

    Optimization.SetPoint = get(h.txtSP, 'String');
    Optimization.ForceStrict = get(h.chkStrict, 'Value');
    Optimization.GenerateReport = get(h.chkSimulate, 'Value');
    Optimization.EnableIterationLimit = get(h.chkLimitIter, 'Value'); % Callback
    Optimization.IterationLimit = get(h.txtLimitIter, 'String');

    % FOMCON specific options
    FOMCON = struct;
    conf   = fomcon('config');
    FOMCON.Time_domain_computations.Performance_index_weight = ...
      conf.FPID_Optimizer.Time_domain_computations.Performance_index_weight;
    FOMCON.Frequency_domain_computations.Num_points_phase = ...
      conf.FPID_Optimizer.Frequency_domain_computations.Num_points_phase;
    FOMCON.Frequency_domain_computations.Phase_comp_weight = ...
      conf.FPID_Optimizer.Frequency_domain_computations.Phase_comp_weight;
    FOMCON.Frequency_domain_computations.W_cg_comp_weight = ...
      conf.FPID_Optimizer.Frequency_domain_computations.W_cg_comp_weight;
    
    % Final structure
    FPID_Config = struct;
    FPID_Config.PlantModel = PlantModel;
    FPID_Config.FPIDParams = FPIDParams;
    FPID_Config.Simulation = Simulation;
    FPID_Config.Optimization = Optimization;
    FPID_Config.FOMCON = FOMCON;

    % Save the structure
    config_struct.FPID_Optimizer_GUI_config = FPID_Config;

end


% --- Executes during object deletion, before destroying properties.
function figFpidOptimTool_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figFpidOptimTool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Restore warning states, if necessary
opt = get(handles.output, 'UserData');
[ex, warns] = cfieldexists(opt, 'warnings_state'); %#ok<ASGLU>
if ~isempty(warns), warning(warns); end


% --- Executes on button press in chkSimulateOnly.
function chkSimulateOnly_Callback(hObject, eventdata, handles)
% hObject    handle to chkSimulateOnly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkSimulateOnly
if get(handles.chkSimulateOnly, 'Value')
    set(handles.chkSimulate, 'Enable', 'Off');
    set(handles.btnOptimize, 'String', 'Simulate');
else
    set(handles.chkSimulate, 'Enable', 'On');
    set(handles.btnOptimize, 'String', 'Optimize');
end

% !-- This is a container function for fpid_optimize used to choose between
% normal optimization and "Simulate only" options
function [Kp, Ki, Kd, lam, mu, Z] = fpid_optim_fpid_optimize(fsim, fopt, G, usesim, handles)

% Check the handles to see if "Simulation only" mode is requested
if get(handles.chkSimulateOnly, 'Value')
    % Get the values to populate the Z structure
    Kp.val  = fopt.p(1);
	Ki.val  = fopt.p(2);
	Kd.val  = fopt.p(3);
	lam.val = fopt.p(4);
	mu.val  = fopt.p(5);

    mindt      = fopt.simtime(1);
    maxtime    = fopt.simtime(3);
    t          = 0:mindt:maxtime;
    
    Z.initialParams = [Kp.val Ki.val Kd.val lam.val mu.val];
    Z.timeSpan      = t;
    
    % Return the same FOPID parameters
    Kp  = Kp.val;
    Ki  = Ki.val;
    Kd  = Kd.val;
    lam = lam.val;
    mu  = mu.val;
    
else
    % Do the optimization
    [Kp, Ki, Kd, lam, mu, Z] = fpid_optimize_(fsim, fopt, G, usesim);
end


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuRefreshSimulinkModelList_Callback(hObject, eventdata, handles)
% hObject    handle to menuRefreshSimulinkModelList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
refreshSimulinkModels(handles);
