function id1 = iddata(id0)
%IDDATA Convert FIDATA object to IDDATA object
% Identification toolbox is required!

% Get the data
id1 = iddata(id0.y, id0.u, id0.dt);

end

