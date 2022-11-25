function varargout = poly2str(varargin)
%POLY2STR Compatibility for older versions of FOMCON; replaced by FPOLY2STR
[varargout{1:nargout}] = fpoly2str(varargin{:});
end

