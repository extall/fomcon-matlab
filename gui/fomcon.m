function varargout = fomcon(sw, config_structure)
%FOMCON Launches the main set of GUIs and governs global configuration
% Usage: fomcon                    launches the main GUI
%        fomcon('config')          launches the toolbox configuration GUI
%        fomcon('config', config)  imports configuration stored in CONFIG
%        config = fomcon('config') returns current configuration structure

    % Define symbolic name for the configuration
    % parameters structure which is kept in the workspace
    config_name = 'fomcon_config';

    % Check input arguments
    if nargin < 1
        % Fallback to main module GUI
        fotf_gui();
    else

        % Only one possible switch now: get/set configuration
        if strcmpi(sw,'config')
            
            % Configuration structure provided, save it
            if nargin == 2 && isstruct(config_structure)
                config = config_structure;
                config = check_compat(config);
                assignin('base', config_name, config_structure);
            else
                % Get the configuration parameters structure
                if varexists(config_name)
                    config = evalin('base', config_name);
                    
                    % Compatibility: check the structure, if some new
                    % fields are missing, add them
                    config = check_compat(config);
                else
                    config = default_config();
                    config = orderfields(config);
                    assignin('base', config_name, config);
                end
                
            end

            % If editing is requested
            if nargout < 1 && ~(nargin == 2)
                new_config = edit_config(config);
                if ~isempty(new_config)
                    assignin('base', config_name, new_config);
                end
            else
                % Otherwise just return the configuration structure
                varargout{1} = config;
            end

        else
            % Wrong switch: launch the main GUI
            fotf_gui();
        end
    end
end


% This function invokes the propertiesGUI() editor
% to change various FOMCON configuration parameters
function new_config = edit_config(config)

    [h, new_config] = propertiesGUI([], config);

end

% For compatibility reasons, check the options structure
function config = check_compat(config)

    % Traverse every field and check if it also exists in the
    % default config. If not, add the default value
    default_config_vals = default_config();
    fields = getfields(default_config_vals);

    for k=1:length(fields)
        if ~cfieldexists(config, fields{k})
            warning(['Found a configuration option that does not exist in present config, setting to default: ', implode(fields{k}, '.')]);
            config = setfield(config, fields{k}{:}, getfield(default_config_vals, fields{k}{:}));
        end
    end
    
    % Order the fields
    config = orderfields(config);

end

% This function returns default toolbox
% configuration parameters for every module
function config = default_config()
    
    % Configuration has several categories
    
    % Core components
    core = struct;
    
    core.General.Model_significant_digits        = 5;
    core.General.Internal_computation_accuracy   = eps;
    
    core.Commensurate_order.Min_comm_order       = 0.01;
    
    core.Frequency_domain.Default_min_freq_exp   = -5;    % Freq. domain
    core.Frequency_domain.Default_max_freq_exp   = +5;
    core.Frequency_domain.Default_num_points     = 1000;

    core.Frequency_domain_computations.Min_freq_exp   = -6;
    core.Frequency_domain_computations.Max_freq_exp   = +6;
    core.Frequency_domain_computations.Num_points     = 500;
    
    core.Root_locus.X_min              = -10;
    core.Root_locus.X_max              = +10;
    core.Root_locus.Y_max              = 10;
    core.Root_locus.Coarse_grid_points = 400;
    core.Root_locus.Fine_grid_points   = 40;
    
    % Approximations
    approximations = struct;
    approximations.Oustaloup_filter.Default_filter_type     = 'oust';
    approximations.Oustaloup_filter.Default_low_freq_bound  = 0.001;
    approximations.Oustaloup_filter.Default_high_freq_bound = 1000;
    approximations.Oustaloup_filter.Default_approx_order    = 5;
    
    % Fractional-order PID controller optimizer
    fpid_optim = struct;
    fpid_optim.Time_domain_computations.Performance_index_weight = 1;
    fpid_optim.Frequency_domain_computations.Num_points_phase    = 250;
    fpid_optim.Frequency_domain_computations.Phase_comp_weight   = 100;
    fpid_optim.Frequency_domain_computations.W_cg_comp_weight    = 100;
    %fpid_optim.Frequency_domain_computations.Gain_margin_limit   = 1000;
    %fpid_optim.Frequency_domain_computations.Phase_margin_limit  = 1000;
    
    % Populate the configuration structure
    config = struct('Core', core, ...
                    'Approximations', approximations, ...
                    'FPID_Optimizer', fpid_optim);
    
end