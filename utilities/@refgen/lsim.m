function varargout = lsim(ref,g)
%LSIM Quick linear simulation using a refgen signal as first argument
%
% Calling sequence: [y, ...] = lsim(ref, g)

varargout{:} = lsim(g,ref.u,ref.t);

end

