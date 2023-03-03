function G = sample(g, method)
%SAMPLE Create a sample member using the specified sampling method

if nargin < 2
    method = 'uniform';
end

% Sample based on method
switch lower(method)
    case 'uniform'
        % Sample values from uncertainty values using uniform dist.
        pb   = g.b;
        pa   = g.a;
        pdel = g.ioDelay;

        b  = pb.a(:,1) + (pb.a(:,2)-pb.a(:,1)).*rand(size(pb.a,1),1);
        nb = pb.na(:,1) + (pb.na(:,2)-pb.na(:,1)).*rand(size(pb.na,1),1);

        a  = pa.a(:,1) + (pa.a(:,2)-pa.a(:,1)).*rand(size(pa.a,1),1);
        na = pa.na(:,1) + (pa.na(:,2)-pa.na(:,1)).*rand(size(pa.na,1),1);

        % Delay
        iodel = [];
        if ~isempty(pdel)
            iodel = pdel(1,1) + rand(1,1) * (pdel(1,2) - pdel(1,1));
        end

        % Create the FOTF object
        G = fotf(a.', na.', b.', nb.', iodel);
        
    otherwise
        error('Unknown sampling method provided.');
end

end
