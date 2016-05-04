function id_t = trim(id,t1,t2)
%TRIM Trims identification dataset from and to the selected time instances
%   ID_T = TRIM(ID, T1, T2) will return the identification dataset ID with
%   from points obtained [T1, T2] and return the resulting set.

    % Check input argument number
    if nargin < 2
       error('FIDATA:TRIM:NotEnoughInputArguments', ...
             'Not enough input arguments.'); 
    end
    
    % Get parameters from identification structure
    t = id.t; dt = id.dt; u = id.u; y = id.y;
    
    % Check condition t1 < t2
    if nargin < 3
        t2 = t(end);
    end
    
    % Time instances cannot be negative
    if t1 < 0 || t2 < 0
        error('FIDATA:TRIM:BadInOutPoints', ...
              'Trim range time instances cannot be negative.');
    end
    
    % T1 must be less or equal to T2
    if t1 > t2
        warning('FIDATA:TRIM:T1GreaterThanT2', ...
                'Incorrect time range: T2 must be greater than T1. Swapping.');
    end
    
    % Find the correct indicies
    t1_ind = find(t>=t1, 1);
    t2_ind = find(t>=t2, 1);
    
    % Trim data
    y_t = y(t1_ind:t2_ind);
    u_t = u(t1_ind:t2_ind);
    t_t = dt*(0:length(y_t)-1);
    
    % Return dataset
    id_t = fidata(y_t, u_t, t_t);
    
end

