function [R,C,L,K,results] = frac_imp_rc_optimize_ph1(params)
%FRAC_IMP_RC_OPTIMIZE_PH1 Obtain R and C component values for an analog
%  approximation of a given fractional-order system

    % Get optimization parameters
    str_fun  = params.structure;     % Used ciruit structure
    orig_fun = params.model;             % Original fractional-order system
    R_series = params.Rser;          % R preferred series
    C_series = params.Cser;          % C preferred series
    w        = params.w;             % Fequency response points
    Nel      = params.Nel;             % Number of discrete components
    Nmc      = params.Nmc;            % Number of Monte-Carlo iterations
    
    % Results array
    results = {};
    
    % No inductors
    L = [];
    
    % Optimization parameters
    % op = optimset('display', 'iter');
    op = set_options('display', 'on');
    
    % Options
    R_MIN = 1;          % 1  Ohm min resistance
    R_MAX = 1e9;        % 1 MOhm max resistance
    C_MIN = 1e-12;      % 1   pF min capacitance
    C_MAX = 1e-3;       % 1   uF max capacitance
    
    % Note here, that resistances and capacitances need to be normed,
    % since during optimization changes need to be more or less of
    % equivalent magnitude. Therefore, we apply the following:
    %
    % R = R / R_MAX
    % C = C / C_MAX
    % thereby bringing exponents closer together.
    
    % Initialize the network, matrix has the form [R1 C1; R2 C2; ...]
    curr_net = [];
    
    % Assign terminating elements
    curr_net(1,1) = params.R1;
    curr_net(1,2) = params.C1;
    
    % Get the ideal response within at the given frequencies
    r_orig = squeeze(freqresp(orig_fun, w));
    
    % Begin optimization
    for k=2:Nel
        % Debug
        disp(['Iteration ' num2str(k-1)]);
        
%         x = optimize(@(x) obj_func(x, str_fun, r_orig, w, curr_net, [R_MAX C_MAX]), ...
%             [0.5 0.5], ...
%             [R_MIN/R_MAX C_MIN/C_MAX], ...
%             [1 1], ...
%             [], [], [], [], [], [], ...
%             op ...
%             );
        
%         x = lsqnonlin(@(x) obj_func(x, str_fun, r_orig, w, curr_net, [R_MAX C_MAX]), ...
%             [0.5 0.5], ...
%             [R_MIN/R_MAX C_MIN/C_MAX], ...
%             [1 1], ...
%             op ...
%             );

        x = GODLIKE(@(x) obj_func(x, str_fun, r_orig, w, curr_net, [R_MAX C_MAX]), ...
            20, ...
            [R_MIN/R_MAX C_MIN/C_MAX], ...
            [1 1], ...
            [], ...
            op ...
            );
        
        % Get the result and save it
        [Rn, Cn] = denormalize_rc(x, [R_MAX C_MAX]);
        
        % Now, convert to preferred series values
        curr_net(k,1) = prefnumerize(Rn,R_series);
        curr_net(k,2) = prefnumerize(Cn,C_series);
        
    end
    
    % Get structure function
    struct_fun = str2func(str_fun);
    
    % Assign R and C arrays
    R = curr_net(:,1);
    C = curr_net(:,2);
    
    % Get the correction gain
    Z = struct_fun(R,C,[],1);
               
    w_av = sqrt(w(1)*w(end));
    K = abs(squeeze(freqresp(orig_fun, w_av))) / ...
        abs(squeeze(freqresp(Z, w_av)));
    
end

function [R,C] = denormalize_rc(x,maxes)
    R  = x(1)*maxes(1);
    C  = x(2)*maxes(2);
end

function f = obj_func(x, str_fun, r_orig, w, curr_net, maxes)

    % Get structure function
    struct_fun = str2func(str_fun);
    
    % Build the RC network
    [R,C] = denormalize_rc(x, maxes);      % Get R,C network values
    
    % Return the impedance
    Z = struct_fun([curr_net(:,1); R], ...
                   [curr_net(:,2); C], ...,
                   [], 1);

    % Constant phase: does not depend on the magnitude
    r_this = squeeze(freqresp(Z, w));
    
    % DEBUG
    if min(size(angle(r_this)) ~= ...
           size(angle(r_orig)))
        r_orig = r_orig';
    end
     
    % Cost function: MSE of phase response
    err = rad2deg(angle(r_orig) - angle(r_this));
    f = sum(err.^2);
    % f = err;

end