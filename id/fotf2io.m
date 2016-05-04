function [K, L, T, y, time] = fotf2io(sim, model, initparams, op)
%FOTF2IO Identify fractional-order tranfer function as integer-order model
%
% Usage: [K,L,T] = FOTF2IO(FSPARAM, MODEL, INITPARAMS, OP)
% where
%        [K, L, T, Y, TIME] - model parameters [K, L, T] 
%                             and original simulation result [Y, T]
%
%        FSPARAM - structure with the following fields (simulation parameters):
%              .plant - FOTF plant to identify
%              .approx - approximation method, should be 'oust' or 'ref'
%              .w - vector with approx. frequencies [wh; wb]
%              .N - approximation order
%
%        MODEL - type of model to identify: 'fopdt', 'ipdt' or 'foipdt'
%
%        INITPARAMS (optional) - initial parameters vector with [K, L, T].
%                                If not provided, defaults to [50, 50, 50]
%
%        OP (optional) - additional optimization parameters for lsqnonlin
%                        defined by OPTIMSET
%
%        See also: fsparam, optimset, lsqnonlin

    % Simulation parameters
    G = sim.plant;
    method = sim.approx;
    wRange = sim.w;
    wb = wRange(1);
    wh = wRange(2);
    N = sim.N;

    if nargin < 4
    
        % Default optimization options
        op = optimset('Display', 'iter');
        
    end
    
    if nargin == 3

        % Initial guess parameters
        K = initparams(1);
        L = initparams(2);
        T = initparams(3);

    elseif nargin < 3 || (nargin > 3 && isempty(initparams))

        % Use defaults
        K = 50;
        L = 50;
        T = 50;

    end
    
    % Get plant approximation
    plant = oustapp(G, wb, wh, N, method);
        
    % Get step response
    [y, t] = step(plant);
    time = t;
    
    % Identify model
    switch lower(model)
       
        case 'fopdt'
            
            x_opt = lsqnonlin(@(x) fopdtfun(x,y,t),...
                                  [K L T], ...
                                  [-Inf 0 0], ...
                                  [Inf Inf Inf], ...
                                  op);
            % Update T parameter
            T = x_opt(3);
        
        case 'ipdt'
            
            x_opt = lsqnonlin(@(x) ipdtfun(x,y,t), ...
                                 [K L], ...
                                 [-Inf 0], ...
                                 [Inf Inf], ...
                                 op);
                             
        case 'foipdt'
            
            x_opt = lsqnonlin(@(x) foipdtfun(x,y,t), ...
                                 [K L T], ...
                                 [-Inf 0 0], ...
                                 [Inf Inf Inf], ...
                                 op);
            T = x_opt(3);
            
        otherwise
            
            error('Unknown model type specified!');
                
    end
    
    % Update parameters
	K = x_opt(1);
	L = x_opt(2);

end

