function varargout = str2poly(varargin)
%STR2POLY Compatibility for older versions of FOMCON; replaced by STR2FPOLY
[varargout{1:nargout}] = str2fpoly(varargin{:});
end

