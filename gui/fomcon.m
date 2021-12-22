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
                assignin('base', config_name, config_structure);
            else
                % Get the configuration parameters structure
                if varexists(config_name)
                    config = evalin('base', config_name);
                else
                    config = default_config();
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


% This function returns default toolbox
% configuration parameters for every module
function config = default_config()
    
    % Configuration has several categories
    
    % Core components
    core = struct;
    
    core.General.Model_significant_digits        = 5;
    core.General.Internal_computation_accuracy   = eps;
    
    core.Commensurate_order.Min_comm_order       = 0.01;
    
    % If high-precision simulation is needed, then change this to 2 or 3
    core.Simulations.GL_Order                    = 1;
    
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