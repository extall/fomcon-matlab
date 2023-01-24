function p1 = sample(p, method)
%SAMPLE Create a sample member using the specified sampling method

if nargin < 2
    method = 'uniform';
end

% Sample based on method
switch lower(method)
    case 'uniform'
        % Sample values from uncertainty values using uniform dist.        
        as  = p.a(:,1) + (p.a(:,2)-p.a(:,1)).*rand(size(p.a,1),1);
        nas = p.na(:,1) + (p.na(:,2)-p.na(:,1)).*rand(size(p.na,1),1);
        p1 = ufpoly([as as], [nas nas], p.symb);
        
    otherwise
        error('Unknown sampling method provided.');
end

end

