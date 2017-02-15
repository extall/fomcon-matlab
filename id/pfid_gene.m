function expr = pfid_gene(nb,na,alpha)
%PFID_GENE Generate a symbolic expression for a candidate FOTF 
%
% This function will create a candidate FO transfer function expression
% for using in parameteric identification using the PFID function.
%
% Calling sequence: PFID_GENE(NB,NA,ALPHA)
% NB - pseudo-order of fractional zero polynomial,
% NA - pseudo-order of fractional pole polynomial,
% ALPHA - commensurate order. If supplied, parameters are not generated
%         for the exponents. Only "p" parameters are populated, while
%         the exponents are generated based on the commensurate order.

    % TODO: make function even more flexible by adding various options
    % for parameter generation, e.g.: provide numerical arrays for
    % coefficients or exponents; allow to use constant parameters;
    % add support for delays.

    % Check input arguments
    if nargin<2
        error('PFID_GENE:NotEnoughInputArguments', ...
            'Not enough input arguments.');
    end
    
    if nargin<3
        alpha = [];
    end

    % Parameter indices
    pz = 1:nb; pp = (nb+1):(nb+na);
    nPz = numel(pz); nPp = numel(pp);
    
    exp_z = ''; exp_p = '';
    
    % Generate the expression
    for k=1:nPz
        exp_z = [exp_z get_fac(pz(k),k,nPz) get_exp(pz(k), alpha, k, nPz)];
    end
    
    for k=1:nPp
        exp_p = [exp_p get_fac(pp(k),k,nPp) get_exp(pp(k), alpha, k, nPp)];
    end
    
    % Final expression
    expr = ['(' exp_z ')/(' exp_p ')'];
    
end

function expr = get_fac(pind, k, n)
    
    if(k==n)
        expr = ['+p' num2str(pind)];
    else
        if (k==1)
            expr = ['p' num2str(pind) '*s'];
        else
            expr = ['+p' num2str(pind) '*s'];
        end
    end
end

function expr = get_exp(pind, alpha, k, n)    
    if (k==n)
        expr = '';
    else
        if isempty(alpha)
            expr = ['^q' num2str(pind)];
        else
            expr = ['^' num2str((n-k)*alpha)];
        end
    end
end

