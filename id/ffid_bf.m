function [q,n,m,G,J] = ffid_bf(idd, method, init, nord, maxiter)
%FFID_BF Best fractional-order model fit for given frequency response.
%
%        Searches for a best-fit fractional model satisfying the frequency-
%        domain data given (extension to FFID). Uses optimize to determine
%        the best parameters
%
% Usage: [Q,N,M,G,J] = FFID_BF(IDD, METHOD, INIT, NORD, MAXITER)
%
% where          Q      - best commensurate order
%                N,M    - pole and zero polynomial orders
%                J      - final error index
%                G      - FOTF object
%
%                IDD        - observed frequency-domain data
%                METHOD     - method used: can be 'l' for Levy or 'v' for
%                             Vinagre
%                INIT       - initial guess in form [Q, N, M]
%                NORD       - highest polynomial order (default: 10)
%                MAXITER    - maximum number of iterations (default:
%                             unlimited)
%
% See also: ffidata, ffid, optimize
    
    % Scaling for pseudo-integer optimization
    scaleFactor = 1e6;

    % Options
    op = optimset('Display', 'iter');
    
    % Check input arguments
    if nargin < 1
        error('FFIDBF:NotEnoughInputArguments', 'Not enough input arguments.');
    end
    
    % Optimization options
    if nargin == 5
        op.MaxIter = maxiter;
    end
    
    % Max order
    if nargin < 4
        nord = 10;
    end
    
    % Initial guess
    if nargin < 3
        iord = 1/scaleFactor;
        init = [1 iord iord];
    end
    
    if nargin < 2
        method = 'l';
    end
    
    % Check method
    if ~(strcmpi(method, 'l') || (strcmpi(method, 'v')))
        warning('FFIDBF:BadMethod', 'Method must be either ''l'' or ''v'': using ''l'' as default.');
        method = 'l';
    end
    
    % Check init
    if isempty(init)
        iord = 1/scaleFactor;
        init = [0.01 iord iord];
    end
    
    % Check nord
    if isempty(nord)
       nord = 10; 
    end
    
    % Scale nord
    nord = nord / scaleFactor;
    
    % Turn off badly scaled matrix warning
    s = warning('off', 'MATLAB:nearlySingularMatrix');

    x_opt = optimize(@(x) ffid_bf_j(x, idd, method, scaleFactor), ...
                      init, ...
                      [0.01 0 0], ...
                      [2 nord nord], ...
                      [], [], [], [], [], [], ...
                      op);
    
    % Return results                  
    [q, n, m] = scale_qnm(x_opt, scaleFactor);
    
    
    if nargout > 3
        [a, na, b, nb, G, J] = ffid(idd, q, [n m], method);
    end
    
    % Restore warning state
    warning(s);

end

% Objective function
function J = ffid_bf_j(x, fd, method, scaleFactor)
    
    [q,n,m] = scale_qnm(x, scaleFactor);
    [a, na, b, nb, G, J] = ffid(fd, q, [n m], method);

end

% Scaling function
function [q,n,m] = scale_qnm(x, scaleFactor)

    q = x(1);
    n = round(x(2)*scaleFactor);
    m = round(x(3)*scaleFactor);
    
end