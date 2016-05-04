function options = set_options(varargin)
% SET_OPTIONS                 Set options for the various optimizers
%
% Usage:
%
%   options = set_options('option1', value1, 'option2', value2, ...)
%
%
%   SET_OPTIONS is an easy way to set all options for the global optimization
%   algorithms PSO, DE, GA, ASA in GODLIKE. All options, and their possible
%   values are:
%
%   ======================================================================
%   General Settings:
%   ======================================================================
%       Display : string, either 'off' (default), 'on' or 'CommandWindow',
%                 'Plot'. This option determines the type of display that
%                 is used to show the algorithm's progress. 'CommandWindow'
%                 (or simply 'on') will show relevant information in the
%                 command window, whereas 'Plot' will make a plot in every
%                 iteration of the current population. Note that 'Plot'
%                 will only work if the number of decision variables is 1
%                 or 2 in case of single-pbjective optimization, or between
%                 1 and 3 objectives for multi-objective optimization.
%                 Please note that using any other display setting than
%                 'off' can significantly slow down the optimization.
%   MaxFunEvals : positive scalar, defining the maximum number of
%                 allowable function evaluations. The default is 100,000.
%                 Note that every objective and constraint function
%                 evaluation will be counted as 1 function evaluation. For
%                 multi-objective optimization, each objective function
%                 will be counted.
%      MaxIters : positive scalar, defining the maximum number of
%                 iterations that can be performed. The default is 20.
%      MinIters : positive scalar. This option defines the minimum amount
%                 of iterations GODLIKE will perform. This is particularly
%                 useful in multi-objective problems with small population
%                 sizes, because this combination increases the probability 
%                 that GODLIKE reports convergence (all fronts are Pareto 
%                 fronts), while a Pareto front of much better quality is 
%                 obtained if some additional shuffles are performed. The
%                 default value is 2. 
%
%   ======================================================================
%   Options specific to the GODLIKE Algorithm:
%   ======================================================================
%       ItersLb : positive scalar. This sets the minimum number of
%                 iterations that will be spent in one of the selected
%                 heuristic optimizers, per GODLIKE iteration. The default
%                 value is 10.
%       ItersUb : positive scalar. This sets the maximum TOTAL amount of
%                 iterations that will be spent in all of the selected
%                 heuristic optimizers combined. The default value is 100.
%
%   ======================================================================
%   General Settings for Single-Objective Optimization:
%   ======================================================================  
%        TolIters: positive scalar. This option defines how many consecutive 
%                  iterations the convergence criteria must hold for each 
%                  individual algorithm, before that algorithm is said to 
%                  have converged. The default setting is 15 iterations. 
%           TolX : positive scalar. Convergence is assumed to be attained, 
%                  if the coordinate differences in all dimensions for a
%                  given amount of consecutive iterations is less than 
%                  [TolX]. This amount of iterations is [TolIters] for each 
%                  individual algorithm, and simply 2 for GODLIKE-iterations. 
%                  The default value is 1e-4.
%         TolFun : positive scalar. Convergence is said to have been 
%                  attained if the value of the objective function decreases 
%                  less than [TolFun] for a given amount of consecutive
%                  iterations. This amount of iterations is [TolIters] for 
%                  each individual algorithm, and simply 2 for the 
%                  GODLIKE-iterations. The default value is 1e-4.
%  AchieveFunVal : scalar. This value is used in conjunction with the
%                  [TolX] and [TolFun] settings. If set, the algorithm will 
%                  FIRST try to achieve this function value, BEFORE enabling
%                  the [TolX] and [TolFun] convergence criteria. By default, 
%                  it is switched off (equal to AchieveFunVal = inf).
%
%   ======================================================================
%   General Settings for Multi-Objective Optimization:
%   ======================================================================
%        SkipTest : If set to 'on', some initial tests that are performed on
%                   the objective and constraint functions. These tests
%                   automatically determine whether the function accepts
%                   vectorized input or not, and how many objectives the
%                   problem has. The default is 'on', but it may be switched
%                   'off'. In case it's switched 'off', the algorithm assumes
%                   all functions accept vectorized input, AND the number of
%                   objectives (the next option) has been given, AND the
%                   dimensionality of the problem is also given (two options
%                   down). The 'off'-switch will be ignored if either of these
%                   demands is not true.
%   NumObjectives : Positive scalar. Sets the number of objectives manually.
%                   When the objective function is a single function that
%                   returns multiple objectives, the algorithm has to first
%                   determine how many objectives there are. This takes some
%                   function evaluations, which may be skipped by setting this
%                   value manually.
%
%   ======================================================================
%   Options specific to the Differential Evolution algorithm:
%   ======================================================================
%          Flb : scalar. This value defines the lower bound for the range
%                from which the scaling parameter will be taken. The
%                default value is -1.5.
%          Fub : scalar. This value defines the upper bound for the range
%                from which the scaling parameter will be taken. The
%                default value is +1.5. These two values may be set equal
%                to each other, in which case the scaling parameter F is
%                simply a constant.
%   CrossConst : positive scalar. It defines the probability with which a
%                new trial individual will be inserted in the new
%                population. The default value is 0.95.
%
%   ======================================================================
%   Options specific to the Genetic Algorithm:
%   ======================================================================
%      Crossprob : positive scalar, defining the probability for crossover
%                  for each individual. The default value is 0.25.
%   MutationProb : positive scalar, defining the mutation probability for
%                  each individual. The default value is 0.1.
%         Coding : string, can either be 'binary' or 'real'. This decides
%                  the coding, or representation, of the variables used by
%                  the genetic algorithm. The default is 'Binary'.
%        NumBits : positive scalar. This options sets the number of bits
%                  to use per decision variable, if the 'Coding' option is
%                  set to 'Binary'. Note that this option is ignored when
%                  the 'Coding' setting is set to 'real'. The default
%                  number of bits is 52 (maximum precision).
%
%   ======================================================================
%   Options specific to the Adaptive Simulated Annealing Algorithm:
%   ======================================================================
%               T0 : positive scalar. This is the initial temperature for 
%                    all particles. If left empty, an optimal one will be
%                    estimated; this is the default. 
%  CoolingSchedule : function handle, with [iteration], [T0], and[T] as
%                    parameters. This function defines the cooling schedule
%                    to be applied each iteration. The default is
%
%                      @(T,T0,iteration) T0 * 0.87^iteration
%
%                    It is only included for completeness, and testing
%                    purposes. Only in rare cases is it beneficial to change
%                    this setting.
%        ReHeating : positive scalar. After an interchange operation in 
%                    GODLIKE, the temperature of an ASA population should
%                    be increased to allow the new individuals to move
%                    over larger portions of the search space. The default
%                    value is 
%
%   ======================================================================
%   Options specific to the Particle Swarm Algorithm:
%   ======================================================================
%           eta1 : scalar < 4. This is the 'social factor', the
%                  acceleration factor in front of the difference with the
%                  particle's position and its neighorhood-best. The
%                  default value is 2. Note that negative values result in
%                  a Repulsive Particle Swarm algorithm.
%           eta2 : scalar < 4. This is the 'cooperative factor', the
%                  acceleration factor in front of the difference with the
%                  particle's position and the location of the global
%                  minimum found so far. The default value is 2.
%           eta3 : scalar < 4. This is the 'nostalgia factor', the
%                  acceleration factor in front of the difference with the
%                  particle's position and its personal-best. The default
%                  value is 0.5.
%          omega : scalar. This is the 'inertial constant', the tendency of
%                  a particle to continue its motion undisturbed. The
%                  default value is 0.5.
%   NumNeighbors : positive scalar. This defines the maximum number of
%                  'neighbors' or 'friends' assigned to each particle. The
%                  default value is 5.
% NetworkTopology: string, equal to either 'fully_connected', 'star', or 
%                  'ring'. This defines the topology of the social network
%                  for each particle. In case 'star' is selected (the 
%                  default), the setting for NumNeighbors will define the 
%                  total number of partiles per star; the same holds in 
%                  case 'ring' is selected. When 'fully_connected' is 
%                  selected however, the value for NumNeighbors will be 
%                  ignored (all particles are connected to all other 
%                  particles). 
%
% see also GODLIKE, pop_multi, pop_single.
    
    % Last edited 30/Jul/2009
    
    % create structure with default values if no input is given
    if (nargin == 0)
        
        % initialize
        options = struct;
        
        % general options
        options.display       = 'off';
        options.MaxFunEvals   = 1e5;
        options.MaxIters      = 20;
        options.MinIters      = 2;        
        options.TolIters      = 15;
        options.TolX          = 1e-4;
        options.TolFun        = 1e-4;
        options.AchieveFunVal = inf;
        
        % function evaluation
        options.num_objectives = 1;
        options.dimensions     = [];
        
        % Differential Evolution
        options.DE.Flb        = -1.5;
        options.DE.Fub        = 1.5;
        options.DE.CrossConst = 0.95;
        
        % genetic algorithm
        options.GA.CrossProb    = 0.5;
        options.GA.MutationProb = 0.1;
        options.GA.Coding       = 'Binary';
        options.GA.NumBits      = 52;
        
        % simulated annealing
        options.ASA.T0           = [];
        options.ASA.CoolingSchedule = @(T, T0, iteration) T0*0.87^iteration;
        options.ASA.ReHeating    = 5;
        
        % particle swarm
        options.PSO.eta1         = 2;
        options.PSO.eta2         = 2;
        options.PSO.eta3         = 0.5;
        options.PSO.omega        = 0.5;
        options.PSO.NumNeighbors = 5;
        options.PSO.NetworkTopology = 'star';
        
        % GODLIKE
        options.GODLIKE.ItersLb = 10;
        options.GODLIKE.ItersUb = 100;
        
        % finished
        return;
        
        % create structure with fields according to user input
    elseif (nargin > 0)
        
        % assign default values
        options = set_options;
        
        % errortrap
        if (mod(nargin, 2) ~= 0)
            error('Please provide values for all the options.')
        end
        
        % loop through all the inputs, and use an "if-else cancer" to
        % create the problem structure
        for i = 1:2:nargin
            option = varargin{i};
            value  = varargin{i+1};
            
            % if option is not recognized, continue to the next argument
            if ~isa(option, 'char')
                throwwarning(option, [], [], []);
                continue;
            end
            
            % General options
            if     strcmpi(option, 'Display')
                if ~ischar(value)
                    throwwarning('Display', 'char', value);
                    continue;
                end
                if     strcmpi(value, 'off')
                    options.display = [];
                elseif strcmpi(value, 'CommandWindow') || strcmpi(value, 'on')
                    options.display = 'CommandWindow';
                elseif strcmpi(value, 'Plot')
                    options.display = 'Plot';
                else
                    error('population:set_options:unknown_display_option',...
                        ['Unknown display type: ', '''', value, '''.'])
                end
            elseif strcmpi(option, 'MaxFunEvals')
                if ~isnumeric(value)
                    throwwarning('MaxFunEvals', 'double', value);
                    continue;
                end
                options.MaxFunEvals = value;
            elseif strcmpi(option, 'MaxIters')
                if ~isnumeric(value)
                    throwwarning('MaxIters', 'double', value);
                    continue;
                end
                options.MaxIters = value;
            elseif strcmpi(option, 'MinIters')
                if ~isnumeric(value)
                    throwwarning('MinIters', 'double', value);
                    continue;
                end
                options.MinIters = value;
            elseif strcmpi(option, 'TolIters')
                if ~isnumeric(value)
                    throwwarning('TolIters', 'double', value);
                    continue;
                end
                options.TolIters = value;
            elseif strcmpi(option, 'TolX')
                if ~isnumeric(value)
                    throwwarning('TolX', 'double', value);
                    continue;
                end
                options.TolX = value;
            elseif strcmpi(option, 'TolFun')
                if ~isnumeric(value)
                    throwwarning('TolFun', 'double', value);
                    continue;
                end
                options.TolFun = value;           
            elseif strcmpi(option, 'AchieveFunVal')
                if ~isnumeric(value)
                    throwwarning('AchieveFunVal', 'double', value);
                    continue;
                end
                options.AchieveFunVal = value;            
                
            % options specific to Differential Evolution
            elseif strcmpi(option, 'Flb')
                if ~isnumeric(value)
                    throwwarning('Flb', 'double', value);
                    continue;
                end
                options.DE.Flb = value;
            elseif strcmpi(option, 'Fub')
                if ~isnumeric(value)
                    throwwarning('Fub', 'double', value);
                    continue;
                end
                options.DE.Fub = value;
            elseif strcmpi(option, 'CrossConst')
                if ~isnumeric(value)
                    throwwarning('CrossConst', 'double', value);
                    continue;
                end
                options.DE.CrossConst = value;
                
            % options specific to Genetic Algorithm
            elseif strcmpi(option, 'MutationProb')
                if ~isnumeric(value)
                    throwwarning('MutationProb', 'double', value);
                    continue;
                end
                options.GA.MutationProb = value;
            elseif strcmpi(option, 'CrossProb')
                if ~isnumeric(value)
                    throwwarning('CrossProb', 'double', value);
                    continue;
                end
                options.GA.CrossProb = value;
            elseif strcmpi(option, 'Coding')
                if ~ischar(value)
                    throwwarning('Coding', 'char', value);
                    continue;
                end
                if     strcmpi(value, 'Real')
                    options.GA.Coding = 'Real';
                elseif strcmpi(value, 'Binary')
                    options.GA.Coding = 'Binary';
                else
                    error('population:set_options:unknown_coding',...
                        ['Unknown coding type: ', '''', value, '''.'])
                end
            elseif strcmpi(option, 'NumBits')
                if ~isnumeric(value)
                    throwwarning('NumBits', 'double', value);
                    continue;
                end
                options.GA.NumBits = value;
                
            % options specific to Adaptive Simulated Annealing
            elseif strcmpi(option, 'T0')
                if ~isnumeric(value)
                    throwwarning('T0', 'double', value);
                    continue;
                end
                options.ASA.T0 = abs(real(value));             
            elseif strcmpi(option, 'CoolingSchedule')
                if ~isa(value, 'function_handle')
                    throwwarning('CoolingSchedule', 'function_handle', value);
                    continue;
                end
                options.ASA.CoolingSchedule = value;
            elseif strcmpi(option, 'ReHeating')
                if ~isa(value, 'double')
                    throwwarning('ReHeating', 'double', value);
                    continue;
                end
                options.ASA.ReHeating = value;
                
            % options specific to Particle Swarm Optimization
            elseif strcmpi(option, 'eta1')
                if ~isnumeric(value)
                    throwwarning('eta1', 'double', value);
                    continue;
                end
                options.PSO.eta1 = value;
            elseif strcmpi(option, 'eta2')
                if ~isnumeric(value)
                    throwwarning('eta2', 'double', value);
                    continue;
                end
                options.PSO.eta2 = value;
            elseif strcmpi(option, 'eta3')
                if ~isnumeric(value)
                    throwwarning('eta3', 'double', value);
                    continue;
                end
                options.PSO.eta3 = value;
            elseif strcmpi(option, 'omega')
                if ~isnumeric(value)
                    throwwarning('omega', 'double', value);
                    continue;
                end
                options.PSO.omega = value;
            elseif strcmpi(option, 'NumNeighbors')
                if ~isnumeric(value)
                    throwwarning('NumNeighbors', 'double', value);
                    continue;
                end
                options.PSO.NumNeighbors = value;
            elseif strcmpi(option, 'NetworkTopology')
                if ~ischar(value)
                    throwwarning('NetworkTopology', 'char', value);
                    continue;
                end
                if strcmpi(value, 'fully_connected')
                    options.PSO.NetworkTopology = 'fully_connected';
                elseif strcmpi(value, 'star')                
                    options.PSO.NetworkTopology = 'star';
                elseif strcmpi(value, 'ring') 
                    options.PSO.NetworkTopology = 'ring';
                else
                    error(['Unknown topology: ''', value, '''.'])
                end
                
            % options specific to GODLIKE algorithm
            elseif strcmpi(option, 'ItersLb')
                if ~isnumeric(value)
                    throwwarning('algiters', 'double', value);
                    continue;
                end
                options.GODLIKE.ItersLb = value;
            elseif strcmpi(option, 'ItersUb')
                if ~isnumeric(value)
                    throwwarning('ItersUb', 'double', value);
                    continue;
                end
                options.GODLIKE.ItersUb = value;
                
                % General Settings
            elseif strcmpi(option, 'SkipTest')
                if ~isnumeric(value)
                    throwwarning('SkipTest', 'char', value);
                    continue;
                end
                if     strcmpi(value, 'on')
                    options.skip_function_test = true;
                elseif strcmpi(value, 'off')
                    options.skip_function_test = false;
                end
            elseif strcmpi(option, 'NumObjectives')
                if ~isnumeric(value)
                    throwwarning('NumObjectives', 'double', value);
                    continue;
                end
                options.num_objectives = value;
            else
                throwwarning(option);
            end % if
        end % for
    end % if
    
    % throw appropriate warning upon abuse of the function
    function throwwarning(option, required, given, varargin)%#ok
        
        % test type
        if nargin == 3
            provided = whos('given');
            provided = provided.class;
            warning('set_options:incorrectvalue', ...
                ['Incorrect class type given for option ''%s'';\n',...
                'required type is ''%s'', received ''%s''.\n',...
                'Using default value...'], option, required, provided);
            
        % unrecognized options will be ignored
        else
            warning('set_options:incorrectoption', ...
                'Unrecognized option, ''%s''. Ignoring it...', num2str(option))
        end % if
    end % nested function
end % function (set options)