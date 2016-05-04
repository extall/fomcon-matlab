function st = imfidata(varargin)
%IMFIDATA Import to FIDATA from another type of data object

% TODO: For now, only support import from a ScopeData structure
% The Scope must have two inputs, the first one collecting OUTPUT data,
% while the second one the INPUT data for the plant of interest.
sd_data = varargin(1);
sd = sd_data{:};

% Regularize the simulation data

% Determine smallest time step
dtmin = min(diff(sd.time));

[y,u,t] = sim_regularize(sd.signals(1).values, ...
                            sd.signals(2).values, ...
                            sd.time, dtmin);
                        
st = fidata(y,u,t);

end

