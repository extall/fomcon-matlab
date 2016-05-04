function params = fpid_recs_ini(varargin)
%FPID_RECS_INI Create a structure with FOPID parameters for FPID_RECS
%
% Usage: PARAMS = FPID_RECS_INI(KP,KI,KD,LAM,MU)
%        PARAMS = FPID_RECS_INI(FPID_VEC)
%
% The output argument PARAMS will then hold the structure with the fields
% as required by the FPID_RECS function
% 
% The input arguments may be supplied by a single vector FPID_VEC
% which has the format [Kp Ki Kd lambda mu].
%
% See also: fpid_recs

params = struct;
if nargin == 1
   vec = varargin{1};
   params.Kp     = vec(1);
   params.Ki     = vec(2);
   params.Kd     = vec(3);
   params.lambda = vec(4);
   params.mu     = vec(5);
elseif nargin == 5
   params.Kp     = varargin{1};
   params.Ki     = varargin{2};
   params.Kd     = varargin{3};
   params.lambda = varargin{4};
   params.mu     = varargin{5}; 
else
    error('FPID_RECS_INI:WrongNumberOfInputArguments', ...
          'Wrong number of input arguments.');
end

end

