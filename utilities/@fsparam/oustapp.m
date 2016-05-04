function Z = oustapp(fsp)
%OUSTAPP Oustaloup approximation for FSPARAM object
%
% Usage: Z = OUSTAPP(FSP)
%
% where  Z is the MATLAB zpk object,
%        FSP is the simulation structure FSPARAM
%
% See also: oustapp, fsparam

    Z = oustapp(fsp.plant, fsp.w(1), fsp.w(2), fsp.N, fsp.approx);

end

