function app = fpid(varargin)
%FPID Fractional-order PID controller design tool.
%
%   This function launches a modern UIFIGURE based interface that replaces
%   the legacy GUIDE implementation. The functional behaviour is preserved
%   while adopting programmatic UI construction to support the latest
%   MATLAB releases.
%
%   app = FPID launches the interface and returns the app controller
%   instance. When an output argument is not requested the app handle is
%   stored internally and the UI lifetime is managed by the figure. The
%   function accepts the same name/value arguments as the historic version:
%   FPID('UserData', SYSNAME) will populate the workspace model field with
%   SYSNAME.
%
%   See also: FOMCON, FRACPID, VAREXISTS, IMPID, IOPID_TUNE.

if nargout
    app = FpidApp(varargin{:});
else
    FpidApp(varargin{:});
end

end

classdef FpidApp < handle
    %FPIDAPP Modern UI controller for the fractional PID design utility.

    properties (Access = private)
        Figure matlab.ui.Figure
        Layout matlab.ui.container.GridLayout
        SystemField matlab.ui.control.EditField
        TimeField matlab.ui.control.EditField
        SvField matlab.ui.control.EditField
        KpField matlab.ui.control.EditField
        KiField matlab.ui.control.EditField
        LambdaField matlab.ui.control.EditField
        KdField matlab.ui.control.EditField
        MuField matlab.ui.control.EditField
        Image matlab.ui.control.Image
        MessageArea matlab.ui.control.TextArea
        ViewButton matlab.ui.control.Button
        SimButton matlab.ui.control.Button
        ExportPlantButton matlab.ui.control.Button
        ExportControllerButton matlab.ui.control.Button
        ApproximateButton matlab.ui.control.Button
        OpenLoopButton matlab.ui.control.Button
        ClosedLoopButton matlab.ui.control.Button
        RealizeButton matlab.ui.control.Button
        ClearButton matlab.ui.control.Button
        ImportButton matlab.ui.control.Button
    end

    methods
        function app = FpidApp(varargin)
            % Constructor builds the UI and applies optional parameters.
            
            % Parse optional legacy parameter/value pairs.
            sysName = '';
            if ~isempty(varargin)
                try
                    sysName = app.parseLegacyArgs(varargin{:});
                catch ME
                    warning('FPID:InvalidArguments', '%s', ME.message);
                end
            end

            buildUI(app);
            populateDefaults(app, sysName);
        end

        function delete(app)
            if ~isempty(app.Figure) && isvalid(app.Figure)
                delete(app.Figure);
            end
        end
    end

    methods (Access = private)
        function sysName = parseLegacyArgs(~, varargin)
            % Supports historic ('UserData', value) invocation style.
            p = inputParser;
            addParameter(p, 'UserData', '', @(x) ischar(x) || isstring(x));
            parse(p, varargin{:});
            sysName = char(p.Results.UserData);
        end

        function buildUI(app)
            app.Figure = uifigure('Name', 'FOMCON - Fractional PID Design', ...
                                  'Position', [100 100 1100 620]);
            app.Figure.CloseRequestFcn = @(src, evt)delete(app);

            app.Layout = uigridlayout(app.Figure, [3 3]);
            app.Layout.RowHeight = {40, '1x', 60};
            app.Layout.ColumnWidth = {280, '1x', 320};
            app.Layout.Padding = [10 10 10 10];
            app.Layout.RowSpacing = 10;
            app.Layout.ColumnSpacing = 12;

            buildHeader(app);
            buildPidPanel(app);
            buildActionPanel(app);
            buildImagePanel(app);
            buildMessageArea(app);
        end

        function buildHeader(app)
            header = uipanel(app.Layout, 'Title', 'Plant / Simulation setup');
            header.Layout.Row = 1;
            header.Layout.Column = [1 3];
            headerGrid = uigridlayout(header, [1 6]);
            headerGrid.ColumnWidth = {120, '1x', 100, 120, 120, 100};
            headerGrid.RowHeight = {'fit'};
            headerGrid.Padding = [10 10 10 10];
            
            uilabel(headerGrid, 'Text', 'Workspace model:', 'HorizontalAlignment', 'right');
            app.SystemField = uieditfield(headerGrid, 'text');
            app.SystemField.ValueChangedFcn = @(s, e)clearMessage(app);

            uilabel(headerGrid, 'Text', 'Time vector:', 'HorizontalAlignment', 'right');
            app.TimeField = uieditfield(headerGrid, 'text');
            app.TimeField.Value = '0:0.1:100';
            app.TimeField.ValueChangedFcn = @(s, e)clearMessage(app);

            uilabel(headerGrid, 'Text', 'Setpoint:', 'HorizontalAlignment', 'right');
            app.SvField = uieditfield(headerGrid, 'text');
            app.SvField.Value = '1';
            app.SvField.ValueChangedFcn = @(s, e)clearMessage(app);
        end

        function buildPidPanel(app)
            pidPanel = uipanel(app.Layout, 'Title', 'Controller parameters');
            pidPanel.Layout.Row = 2;
            pidPanel.Layout.Column = 1;

            grid = uigridlayout(pidPanel, [5 2]);
            grid.RowHeight = repmat({34}, 1, 5);
            grid.ColumnWidth = {80, '1x'};
            grid.Padding = [10 10 10 10];
            grid.RowSpacing = 8;

            uilabel(grid, 'Text', 'Kp:', 'HorizontalAlignment', 'right');
            app.KpField = uieditfield(grid, 'text', 'Value', '1');
            app.KpField.ValueChangedFcn = @(s, e)clearMessage(app);

            uilabel(grid, 'Text', 'Ki:', 'HorizontalAlignment', 'right');
            app.KiField = uieditfield(grid, 'text', 'Value', '1');
            app.KiField.ValueChangedFcn = @(s, e)clearMessage(app);

            uilabel(grid, 'Text', 'λ (Integral order):', 'HorizontalAlignment', 'right');
            app.LambdaField = uieditfield(grid, 'text', 'Value', '0.5');
            app.LambdaField.ValueChangedFcn = @(s, e)clearMessage(app);

            uilabel(grid, 'Text', 'Kd:', 'HorizontalAlignment', 'right');
            app.KdField = uieditfield(grid, 'text', 'Value', '1');
            app.KdField.ValueChangedFcn = @(s, e)clearMessage(app);

            uilabel(grid, 'Text', 'μ (Derivative order):', 'HorizontalAlignment', 'right');
            app.MuField = uieditfield(grid, 'text', 'Value', '0.5');
            app.MuField.ValueChangedFcn = @(s, e)clearMessage(app);
        end

        function buildActionPanel(app)
            actionPanel = uipanel(app.Layout, 'Title', 'Actions');
            actionPanel.Layout.Row = 2;
            actionPanel.Layout.Column = 2;

            grid = uigridlayout(actionPanel, [6 2]);
            grid.RowHeight = repmat({32}, 1, 6);
            grid.ColumnWidth = {'1x', '1x'};
            grid.RowSpacing = 8;
            grid.Padding = [10 10 10 10];

            app.ViewButton = uibutton(grid, 'Text', 'View controller', ...
                'ButtonPushedFcn', @(src, evt)onView(app));
            app.ViewButton.Layout.Row = 1;
            app.ViewButton.Layout.Column = 1;

            app.SimButton = uibutton(grid, 'Text', 'Simulate', ...
                'ButtonPushedFcn', @(src, evt)onSimulate(app));
            app.SimButton.Layout.Row = 1;
            app.SimButton.Layout.Column = 2;

            app.OpenLoopButton = uibutton(grid, 'Text', 'Open-loop Bode', ...
                'ButtonPushedFcn', @(src, evt)onOpenLoopBode(app));
            app.OpenLoopButton.Layout.Row = 2;
            app.OpenLoopButton.Layout.Column = 1;

            app.ClosedLoopButton = uibutton(grid, 'Text', 'Closed-loop Bode', ...
                'ButtonPushedFcn', @(src, evt)onClosedLoopBode(app));
            app.ClosedLoopButton.Layout.Row = 2;
            app.ClosedLoopButton.Layout.Column = 2;

            app.ExportControllerButton = uibutton(grid, 'Text', 'Export control system', ...
                'ButtonPushedFcn', @(src, evt)onExportControlSystem(app));
            app.ExportControllerButton.Layout.Row = 3;
            app.ExportControllerButton.Layout.Column = 1;

            app.ExportPlantButton = uibutton(grid, 'Text', 'Export PID', ...
                'ButtonPushedFcn', @(src, evt)onExportPID(app));
            app.ExportPlantButton.Layout.Row = 3;
            app.ExportPlantButton.Layout.Column = 2;

            app.ApproximateButton = uibutton(grid, 'Text', 'Approximate PID', ...
                'ButtonPushedFcn', @(src, evt)onApproximate(app));
            app.ApproximateButton.Layout.Row = 4;
            app.ApproximateButton.Layout.Column = 1;

            app.RealizeButton = uibutton(grid, 'Text', 'Realize controller', ...
                'ButtonPushedFcn', @(src, evt)onRealize(app));
            app.RealizeButton.Layout.Row = 4;
            app.RealizeButton.Layout.Column = 2;

            app.ImportButton = uibutton(grid, 'Text', 'Import PID config', ...
                'ButtonPushedFcn', @(src, evt)onImportConfig(app));
            app.ImportButton.Layout.Row = 5;
            app.ImportButton.Layout.Column = 1;

            app.ClearButton = uibutton(grid, 'Text', 'Clear fields', ...
                'ButtonPushedFcn', @(src, evt)onClear(app));
            app.ClearButton.Layout.Row = 5;
            app.ClearButton.Layout.Column = 2;

            % spacer row for layout aesthetic
            grid.RowHeight{6} = '1x';
        end

        function buildImagePanel(app)
            imgPanel = uipanel(app.Layout, 'Title', 'Architecture');
            imgPanel.Layout.Row = [1 2];
            imgPanel.Layout.Column = 3;
            imgPanel.Scrollable = 'on';

            imgGrid = uigridlayout(imgPanel, [1 1]);
            imgGrid.Padding = [10 10 10 10];

            app.Image = uiimage(imgGrid);
            imgPath = fullfile(fileparts(mfilename('fullpath')), 'fpid.jpg');
            if exist(imgPath, 'file')
                app.Image.ImageSource = imgPath;
            end
            app.Image.ScaleMethod = 'fit';
        end

        function buildMessageArea(app)
            msgPanel = uipanel(app.Layout, 'Title', 'Messages');
            msgPanel.Layout.Row = 3;
            msgPanel.Layout.Column = [1 3];
            msgPanel.Padding = [10 10 10 10];

            app.MessageArea = uitextarea(msgPanel, 'Editable', 'off');
            app.MessageArea.Value = {'Welcome to the fractional PID design tool.'};
            app.MessageArea.FontName = 'monospaced';
            app.MessageArea.Layout.Row = 1;
            app.MessageArea.Layout.Column = 1;
        end

        function populateDefaults(app, sysName)
            if ~isempty(sysName)
                app.SystemField.Value = sysName;
            end
        end

        function clearMessage(app)
            if isvalid(app.MessageArea)
                app.MessageArea.Value = {''};
            end
        end

        function appendMessage(app, msg)
            if ~isvalid(app.MessageArea)
                return;
            end
            current = app.MessageArea.Value;
            current(end+1) = {char(msg)}; %#ok<AGROW>
            app.MessageArea.Value = current;
            drawnow limitrate;
        end

        function [modelName, plant] = fetchPlant(app)
            modelName = strtrim(app.SystemField.Value);
            if isempty(modelName)
                error('FPID:NoPlant', 'Specify a workspace model name first.');
            end
            [existsFlag, classMatch, cls] = varexists(modelName);
            if ~existsFlag
                error('FPID:MissingPlant', 'Object %s does not exist in the workspace.', modelName);
            end
            if ~classMatch
                error('FPID:InvalidPlant', 'Object %s is not a supported LTI or FOTF model.', modelName);
            end
            plant = evalin('base', modelName);
            appendMessage(app, sprintf('Using plant ''%s'' (%s).', modelName, cls));
        end

        function [Kp, Ki, lambda, Kd, mu] = getPidParams(app)
            [Kp, Ki, lambda, Kd, mu] = evaluateFields(app, ...
                app.KpField, app.KiField, app.LambdaField, app.KdField, app.MuField);
        end

        function [Kp, Ki, lambda, Kd, mu] = evaluateFields(app, varargin)
            vals = cell(1, numel(varargin));
            defaults = {1, 1, 0.5, 1, 0.5};
            for idx = 1:numel(varargin)
                field = varargin{idx};
                txt = strtrim(field.Value);
                if isempty(txt)
                    field.Value = num2str(defaults{idx});
                    vals{idx} = defaults{idx};
                else
                    vals{idx} = evalin('base', txt);
                end
            end
            Kp = vals{1};
            Ki = vals{2};
            lambda = vals{3};
            Kd = vals{4};
            mu = vals{5};
        end

        function [timeVec, sv] = getSimulationSetup(app)
            timeStr = strtrim(app.TimeField.Value);
            if isempty(timeStr)
                timeStr = '0:0.1:30';
                app.TimeField.Value = timeStr;
            end
            timeVec = evalin('base', timeStr);

            svStr = strtrim(app.SvField.Value);
            if isempty(svStr)
                svStr = '1';
                app.SvField.Value = svStr;
            end
            svVal = evalin('base', svStr);
            sv = svVal * ones(size(timeVec));
        end

        function pid = buildPid(app)
            [Kp, Ki, lambda, Kd, mu] = getPidParams(app);
            pid = fracpid(Kp, Ki, lambda, Kd, mu);
        end

        function [pid, pidType] = buildPidWithType(app)
            [Kp, Ki, lambda, Kd, mu] = getPidParams(app);
            [pid, pidType] = fracpid(Kp, Ki, lambda, Kd, mu);
        end

        function onView(app)
            try
                [pid, pidType] = buildPidWithType(app);
                appendMessage(app, sprintf('Current controller (%s):', pidType));
                disp(pid);
            catch ME
                handleError(app, ME);
            end
        end

        function onSimulate(app)
            try
                [modelName, plant] = fetchPlant(app);
                [timeVec, sv] = getSimulationSetup(app);
                [pid, pidType] = buildPidWithType(app);

                controller = pid;
                if isa(plant, 'fotf')
                    delay = plant.ioDelay;
                elseif isprop(plant, 'ioDelay')
                    delay = plant.ioDelay;
                else
                    delay = 0;
                end

                if any(delay(:))
                    controller = oustapp(pid, 1e-4, 1e4, 5, 'oust');
                    appendMessage(app, 'Delay detected. Using Oustaloup approximation for simulation.');
                    if ~isproper(controller)
                        controller = toproper(controller, 1e4);
                    end
                end

                closedLoop = feedback(controller * plant, 1);

                fig = uifigure('Name', sprintf('%s control response', pidType));
                ax = uiaxes(fig);
                y = lsim(closedLoop, sv, timeVec);
                plot(ax, timeVec, y, 'LineWidth', 1.4);
                hold(ax, 'on');
                plot(ax, timeVec, sv, '--r', 'LineWidth', 1);
                hold(ax, 'off');
                grid(ax, 'on');
                xlabel(ax, 'Time [s]');
                ylabel(ax, 'Amplitude');
                legend(ax, {'Response', 'Setpoint'}, 'Location', 'best');
                title(ax, sprintf('%s control system response for %s', pidType, modelName));
            catch ME
                handleError(app, ME);
            end
        end

        function onExportPID(app)
            try
                prompt = {'Workspace variable name:'};
                answer = inputdlg(prompt, 'Export PID-FOTF to Workspace', 1, {''}, struct('WindowStyle', 'normal'));
                if isempty(answer)
                    return;
                end
                varName = strtrim(answer{1});
                if isempty(varName)
                    error('FPID:InvalidName', 'Specify a valid workspace variable name.');
                end
                pid = buildPid(app);
                assignin('base', varName, pid);
                appendMessage(app, sprintf('PID controller exported to workspace variable ''%s''.', varName));
            catch ME
                handleError(app, ME);
            end
        end

        function onExportControlSystem(app)
            try
                [modelName, plant] = fetchPlant(app);
                [pid, pidType] = buildPidWithType(app);

                if isa(plant, 'fotf')
                    controller = pid;
                    prompts = {'Workspace variable name:'};
                    defaults = {sprintf('%s_control', modelName)};
                    answer = inputdlg(prompts, 'Save control system', 1, defaults, ...
                        struct('WindowStyle', 'normal', 'Resize', 'on'));
                    if isempty(answer)
                        return;
                    end
                    varName = strtrim(answer{1});
                    if isempty(varName)
                        varName = defaults{1};
                    end
                    fullCtrl = feedback(controller * plant, 1);
                    assignin('base', varName, fullCtrl);
                    appendMessage(app, sprintf('Control system exported to workspace variable ''%s''.', varName));
                    appendMessage(app, sprintf('Controller type: %s', pidType));
                    return;
                end

                prompts = {'Workspace variable name:', ...
                           'Approximation type (oust or ref):', ...
                           'Low frequency bound wb [rad/s]:', ...
                           'High frequency bound wh [rad/s]:', ...
                           'Order of approximation:'};
                defaults = {sprintf('%s_control', modelName), 'oust', '0.0001', '10000', '5'};
                answer = inputdlg(prompts, 'Export parameters', 1, defaults, ...
                    struct('WindowStyle', 'normal', 'Resize', 'on'));
                if isempty(answer)
                    return;
                end
                ctrlVar = strtrim(answer{1});
                if isempty(ctrlVar)
                    ctrlVar = defaults{1};
                end
                ofType = answer{2};
                ofWb = str2double(answer{3});
                ofWh = str2double(answer{4});
                ofN = str2double(answer{5});
                controller = oustapp(pid, ofWb, ofWh, ofN, ofType);
                if ~isproper(controller)
                    controller = toproper(controller, ofWh);
                end
                plant = ss(plant);
                fullCtrl = feedback(controller * plant, 1);
                assignin('base', ctrlVar, fullCtrl);
                appendMessage(app, sprintf('Approximate control system stored to ''%s''.', ctrlVar));
            catch ME
                handleError(app, ME);
            end
        end

        function onApproximate(app)
            try
                [modelName, plant] = fetchPlant(app);
                [pid, pidType] = buildPidWithType(app);
                paramsPrompt = {'Approximation type (oust or ref):', ...
                                'Low frequency bound wb [rad/s]:', ...
                                'High frequency bound wh [rad/s]:', ...
                                'Order of approximation:'};
                defaults = {'oust', '0.0001', '10000', '5'};
                answer = inputdlg(paramsPrompt, 'PID approximation parameters', 1, defaults, ...
                    struct('WindowStyle', 'normal', 'Resize', 'on'));
                if isempty(answer)
                    return;
                end
                ofType = answer{1};
                ofWb = str2double(answer{2});
                ofWh = str2double(answer{3});
                ofN = str2double(answer{4});
                pidApprox = oustapp(pid, ofWb, ofWh, ofN, ofType);
                if ~isproper(pidApprox)
                    pidApprox = toproper(pidApprox, ofWh);
                end
                assignin('base', sprintf('%s_pidApprox', modelName), pidApprox);
                appendMessage(app, sprintf('Stored %s PID approximation in workspace.', pidType));
            catch ME
                handleError(app, ME);
            end
        end

        function onOpenLoopBode(app)
            try
                [modelName, plant] = fetchPlant(app);
                [pid, pidType] = buildPidWithType(app);
                defaultFreq = 'logspace(-5,5,1000)';
                if isa(plant, 'fotf')
                    freqAnswer = inputdlg({'Frequencies of interest [rad/s]:'}, ...
                        'Bode plot parameters', 1, {defaultFreq}, ...
                        struct('WindowStyle', 'normal', 'Resize', 'on'));
                    if isempty(freqAnswer)
                        return;
                    end
                    w_ev = evalin('base', freqAnswer{1});
                    controller = pid;
                    plantToUse = plant;
                else
                    paramsPrompt = {'Approximation type (oust or ref):', ...
                                    'Low frequency bound wb [rad/s]:', ...
                                    'High frequency bound wh [rad/s]:', ...
                                    'Order of approximation:', ...
                                    'Frequencies of interest [rad/s]:'};
                    defaults = {'oust', '0.0001', '10000', '5', defaultFreq};
                    answer = inputdlg(paramsPrompt, 'Bode plot parameters', 1, defaults, ...
                        struct('WindowStyle', 'normal', 'Resize', 'on'));
                    if isempty(answer)
                        return;
                    end
                    ofType = answer{1};
                    ofWb = str2double(answer{2});
                    ofWh = str2double(answer{3});
                    ofN = str2double(answer{4});
                    w_ev = evalin('base', answer{5});
                    controller = oustapp(pid, ofWb, ofWh, ofN, ofType);
                    if ~isproper(controller)
                        controller = toproper(controller, ofWh);
                    end
                    plantToUse = ss(plant);
                end

                fig = figure('Name', sprintf('%s * %s open-loop response', pidType, modelName));
                bode(controller * plantToUse, w_ev);
                grid on;
            catch ME
                handleError(app, ME);
            end
        end

        function onClosedLoopBode(app)
            try
                [modelName, plant] = fetchPlant(app);
                [pid, pidType] = buildPidWithType(app);
                defaultFreq = 'logspace(-5,5,1000)';
                if isa(plant, 'fotf')
                    answer = inputdlg({'Frequencies of interest [rad/s]:'}, ...
                        'Bode plot parameters', 1, {defaultFreq}, ...
                        struct('WindowStyle', 'normal', 'Resize', 'on'));
                    if isempty(answer)
                        return;
                    end
                    w_ev = evalin('base', answer{1});
                    controller = pid;
                else
                    paramsPrompt = {'Approximation type (oust or ref):', ...
                                    'Low frequency bound wb [rad/s]:', ...
                                    'High frequency bound wh [rad/s]:', ...
                                    'Order of approximation:', ...
                                    'Frequencies of interest [rad/s]:'};
                    defaults = {'oust', '0.0001', '10000', '5', defaultFreq};
                    answer = inputdlg(paramsPrompt, 'Bode plot parameters', 1, defaults, ...
                        struct('WindowStyle', 'normal', 'Resize', 'on'));
                    if isempty(answer)
                        return;
                    end
                    ofType = answer{1};
                    ofWb = str2double(answer{2});
                    ofWh = str2double(answer{3});
                    ofN = str2double(answer{4});
                    w_ev = evalin('base', answer{5});
                    controller = oustapp(pid, ofWb, ofWh, ofN, ofType);
                    if ~isproper(controller)
                        controller = toproper(controller, ofWh);
                    end
                    plant = ss(plant);
                end
                closedLoop = feedback(controller * plant, 1);
                fig = figure('Name', sprintf('%s closed-loop frequency response', modelName));
                bode(closedLoop, w_ev);
                grid on;
            catch ME
                handleError(app, ME);
            end
        end

        function onRealize(app)
            try
                params.Kp = app.KpField.Value;
                params.Ki = app.KiField.Value;
                params.Kd = app.KdField.Value;
                params.lam = app.LambdaField.Value;
                params.mu = app.MuField.Value;
                impid('UserData', params);
            catch ME
                handleError(app, ME);
            end
        end

        function onClear(app)
            app.SystemField.Value = '';
            app.KpField.Value = '1';
            app.KiField.Value = '1';
            app.LambdaField.Value = '0.5';
            app.KdField.Value = '1';
            app.MuField.Value = '0.5';
            app.TimeField.Value = '0:0.1:100';
            app.SvField.Value = '1';
            app.MessageArea.Value = {'Fields reset to defaults.'};
        end

        function onImportConfig(app)
            try
                [filename, path] = uigetfile({'*.mat', 'MAT-files (*.mat)'}, 'Import PID configuration');
                if isequal(filename, 0) || isequal(path, 0)
                    return;
                end
                fileData = load(fullfile(path, filename));
                if ~isfield(fileData, 'FPID_Optimizer_GUI_config')
                    error('FPID:InvalidFile', 'Selected file does not contain FPID configuration.');
                end
                cfg = fileData.FPID_Optimizer_GUI_config;
                app.KpField.Value = cfg.FPIDParams.Kp;
                app.KiField.Value = cfg.FPIDParams.Ki;
                app.KdField.Value = cfg.FPIDParams.Kd;
                app.LambdaField.Value = cfg.FPIDParams.Lam;
                app.MuField.Value = cfg.FPIDParams.Mu;
                appendMessage(app, sprintf('PID parameters imported from %s.', filename));
            catch ME
                handleError(app, ME);
            end
        end

        function handleError(app, ME)
            appendMessage(app, sprintf('Error: %s', ME.message));
            errordlg(ME.message, 'FOMCON FPID error', 'modal');
        end
    end
end
