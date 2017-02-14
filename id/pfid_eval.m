function G = pfid_eval(expr, params)
%PFID_EVAL Evaluate EXPR under PARAMS as returned by PFID and PFID_

% Check input arguments
if nargin < 2
    error('Not enough input arguments.');
end

% Iterate over fields and create the arrays
form = '([p|q])_?([0-9]+)';

% Permissable parameter names are p and q.
% No other parameters are read.
my_fields = fieldnames(params);

% Make p and q vectors at least the same length as parameter number
numfields = numel(my_fields);
p = zeros(numfields,1); q = zeros(numfields,1);

% Now, iterate through all fields and assign correct indices
for k=1:numfields
   
   % Value of parameter 
   nowp = params.(my_fields{k});
    
   % Parse parameter name and index
   tk = regexp(my_fields{k}, form, 'tokens');
   tk_n = tk{1}{1};
   tk_i = str2num(tk{1}{2});
   
   switch tk_n
       case 'p'
           p(tk_i) = nowp;
       case 'q'
           q(tk_i) = nowp;
       otherwise
           warning(['Orphaned parameter ' tk_n 'in parameter structure.']);
   end

end

% Replace id parameters with vectors
expr      = regexprep(expr, form, '$1($2)');

% Finally, evaluate the expression
s = fotf('s');
G = eval(expr);

end

